#include "precomp.h"
#include "mygame.h"


using namespace openvdb::math;






Game* game = new MyGame();

struct _cam
{
	mat4 ct;
#define cam_x ct.cell[3]
#define cam_y ct.cell[7]
#define cam_z ct.cell[11]

	void CalculateMatrix(float3& x, float3& y, float3& z)
	{
		x = make_float3(ct.cell[0], ct.cell[4], ct.cell[8]);
		y = make_float3(ct.cell[1], ct.cell[5], ct.cell[9]);
		z = make_float3(ct.cell[2], ct.cell[6], ct.cell[10]);
	}
	void TranslateRelative(float3 T)
	{
		float3 x, y, z;
		CalculateMatrix(x, y, z);
		float3 delta = T.x * x + T.y * y + T.z * z;
		float3 trans = make_float3(cam_x, cam_y, cam_z) + delta;
		cam_x = trans.x, cam_y = trans.y, cam_z = trans.z;
	}
	void TranslateTarget(float3 T)
	{
		float3 x, y, z;
		CalculateMatrix(x, y, z);
		float3 delta = T.x * x + T.y * y + T.z * z;
		z = normalize(z + delta);
		x = normalize(cross(z, normalize(make_float3(0, 1, 0))));
		y = cross(x, z);
		ct[0] = x.x, ct[4] = x.y, ct[8] = x.z;
		ct[1] = y.x, ct[5] = y.y, ct[9] = y.z;
		ct[2] = z.x, ct[6] = z.y, ct[10] = z.z;
	}


	std::bitset<GLFW_KEY_LAST> keystates;

	void update(float dt)
	{
		float rot = dt * 0.00175f;
		float spd = dt * (keystates[GLFW_KEY_LEFT_SHIFT] ? 1.0f : 0.1f);
		if (keystates[GLFW_KEY_A]) TranslateRelative(make_float3(-spd, 0, 0));
		if (keystates[GLFW_KEY_D]) TranslateRelative(make_float3(spd, 0, 0));
		if (keystates[GLFW_KEY_W]) TranslateRelative(make_float3(0, 0, spd));
		if (keystates[GLFW_KEY_S]) TranslateRelative(make_float3(0, 0, -spd));
		if (keystates[GLFW_KEY_E]) TranslateRelative(make_float3(0, spd, 0));
		if (keystates[GLFW_KEY_Q]) TranslateRelative(make_float3(0, -spd, 0));

		if (keystates[GLFW_KEY_UP]) TranslateTarget(make_float3(0, -rot, 0));
		if (keystates[GLFW_KEY_DOWN]) TranslateTarget(make_float3(0, rot, 0));
		if (keystates[GLFW_KEY_LEFT]) TranslateTarget(make_float3(-rot, 0, 0));
		if (keystates[GLFW_KEY_RIGHT]) TranslateTarget(make_float3(rot, 0, 0));
	}

} cam;

openvdb::FloatGrid::Ptr Load_cached(const char* file)
{
	uint64_t key = (uint64_t)MAPWIDTH << 32 | MAPHEIGHT << 16 | MAPDEPTH;
	
	std::stringstream stream;
	stream << file << "." << key << ".vdb";
	std::string cachedF(stream.str());
	//std::cout << "trying " << cachedF <<std::endl;

	if (FileExists(cachedF.c_str()))
	{
		openvdb::io::File newFile(cachedF);
		newFile.open();
		openvdb::GridBase::Ptr baseGrid = newFile.getGrids()->at(0);
		openvdb::FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);
		newFile.close();

		CoordBBox bbox = grid->evalActiveVoxelBoundingBox();
		//std::cout << "scale" << bbox.dim() << std::endl;

		return grid;
	}
	else
	{
		openvdb::io::File newFile(file);
		newFile.open();
		openvdb::GridBase::Ptr baseGrid = newFile.getGrids()->at(0);
		openvdb::FloatGrid::Ptr gridIN = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);
		newFile.close();

		// rescale into map dims (if needed)
		CoordBBox bbox = gridIN->evalActiveVoxelBoundingBox();
		std::cout << "scale" << bbox.dim() << std::endl;
		float w = bbox.dim().x();
		float h = bbox.dim().y();
		float d = bbox.dim().z();
		w = MAPWIDTH < w ? (float)MAPWIDTH / w : 1;
		h = MAPHEIGHT < h ? (float)MAPHEIGHT / h : 1;
		d = MAPDEPTH < d ? (float)MAPDEPTH / d : 1;
		std::cout << "trying with scale " << w << " " << h << " " << d << " res:" << std::min(std::min(w, h), d) << std::endl;

		openvdb::FloatGrid::Ptr grid;
		float scale = 0.9 * std::min(std::min(w, h), d);
		{
			Timer t;
			grid = openvdb::FloatGrid::create();

			Vec3d trans;
			trans.x() = (MAPWIDTH / 2);
			trans.z() = (MAPDEPTH / 2);
			trans.y() = (MAPHEIGHT/2);

			openvdb::tools::GridTransformer transformer(Vec3s::zero(), Vec3d(scale, scale, scale), Vec3s::zero(), trans);
			transformer.transformGrid<openvdb::tools::QuadraticSampler, openvdb::FloatGrid>(*gridIN, *grid);

			grid->tree().prune();

		//	grid->setTransform(gridIN->transformPtr());

			std::cout << " transform " << t.elapsed() << std::endl;
		}

		
		openvdb::io::File newFile2(cachedF);
		newFile2.write({ grid });
		newFile2.close();
		

		return grid;
	}	
}

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------

float Remap(float value, float from1, float to1, float from2, float to2) {
	return (value - from1) / (to1 - from1) * (to2 - from2) + from2;
}

template<typename GridType>
struct LeafProcessor
{
	using TreeType = typename GridType::TreeType;
	using LeafNode = typename TreeType::ValueOnCIter;
	using IterRange = openvdb::tree::IteratorRange<typename TreeType::ValueOnCIter>;
	void operator()(IterRange& range) const
	{
		World* world = GetWorld();

		// Note: this code must be thread-safe.
		// Iterate over a subrange of the leaf iterator's iteration space.
		for (; range; ++range) {
			// Retrieve the leaf node to which the iterator is pointing.
			const openvdb::FloatGrid::ValueOnCIter& iter = range.iterator();
			// Update the global counter.
			//activeLeafVoxelCount.fetch_and_add(leaf.onVoxelCount());

			auto c = iter.getCoord();
			auto ii = c.asVec3I();
			int x = ii.x();
			int y = ii.y();
			int z = ii.z();
			unsigned char vc = LIGHTBLUE;
			world->Set(x, y, z, vc);
		}
	}
};

void MyGame::Init()
{
	cam.ct = mat4::LookAt(make_float3(MAPWIDTH, MAPHEIGHT *.66f, -MAPDEPTH*.66f), make_float3(MAPWIDTH/2, MAPHEIGHT/2, MAPDEPTH/2));

	// This function is for functionality that you only want to run once,
	// at the start of your game. You can setup the scene here, load some
	// sprites, and so on.

	openvdb::initialize();
	openvdb::FloatGrid::Ptr grid = Load_cached("C:/Code/vox/bunny_cloud.vdb");


#if 0
	World* world = GetWorld();	
	{
		Timer t;
		using FloatLeafProc = LeafProcessor<openvdb::FloatGrid>;
		FloatLeafProc proc;
		FloatLeafProc::IterRange range(grid->tree().cbeginValueOn());
		tbb::parallel_for(range, proc);
		std::cout << " inject MT " << t.elapsed() << std::endl;
	}
#else
	// inject into map
	{
		Timer t;

		World* world = GetWorld();
		int iii = 0;
		for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter; ++iter)
		{
			auto c = iter.getCoord();
			auto ii = c.asVec3I();
			int x = ii.x() ;//+ MAPWIDTH / 2;
			int y = ii.y() ;//+ MAPHEIGHT / 2;
			int z = ii.z() ;//+ MAPDEPTH / 2;
			//	if (x > 0 && y > 0 && z > 0)
			//		if (x < MAPWIDTH  && y < MAPHEIGHT && z < MAPDEPTH)
			{
				float v = *iter;
				//v *= 0.5;
				//v += 0.5;

				//v = Remap(v, -.3f, .3f, 0, 1);

				unsigned char vc = LIGHTBLUE;// (unsigned char)(max(min(255, (int)(v * 256.f)), 1));

				world->Set(x, y, z, vc);
			//	std::cout << c << " " << *iter << std::endl;
			}
		}
		std::cout << " inject ST " << t.elapsed() << std::endl;
	}
#endif

 }



void MyGame::KeyUp(int key)
{
	cam.keystates[key] = 0;
}
void MyGame::KeyDown(int key)
{
	cam.keystates[key] = 1;
}

// -----------------------------------------------------------
// Main application tick function
// -----------------------------------------------------------
void MyGame::Tick( float deltaTime )
{
	cam.update(deltaTime);
	// This function gets called once per frame by the template code.
	World* world = GetWorld();
	

	world->SetCameraMatrix(cam.ct);





}