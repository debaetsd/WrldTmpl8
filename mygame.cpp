#include "precomp.h"
#include "mygame.h"
using namespace openvdb::math;

Game* game = new MyGame();
World* world = GetWorld();

namespace util
{
	static char _filePickerScratch[MAX_PATH];
	const char* Filepicker()
	{
		_filePickerScratch[0] = '\0';

		OPENFILENAME ofn;
		HWND hwnd = glfwGetWin32Window(GetWindow());

		// Initialize OPENFILENAME
		ZeroMemory(&ofn, sizeof(ofn));
		ofn.lStructSize = sizeof(ofn);
		ofn.hwndOwner = hwnd;
		ofn.lpstrFile = _filePickerScratch;
		ofn.nMaxFile = MAX_PATH;
		ofn.lpstrFilter = "openvdb\0*.vdb\0\0";
		ofn.nFilterIndex = 1;
		ofn.lpstrFileTitle = NULL;
		ofn.nMaxFileTitle = 0;
		ofn.lpstrInitialDir = NULL;
		ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

		if (GetOpenFileName(&ofn) == TRUE)
		{
			return _filePickerScratch;
		}
		return NULL;
	}

	// maths stolen from lighthouse 2 (https://github.com/jbikker/lighthouse2)
	struct camera
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
		std::bitset<GLFW_MOUSE_BUTTON_LAST> mouseStates;
		float2 mousepos, lastmouse;

		void update(float dt)
		{
			if (mouseStates[GLFW_MOUSE_BUTTON_RIGHT])
			{
				float spd = dt * (keystates[GLFW_KEY_LEFT_SHIFT] ? 1.0f : 0.1f);
				if (keystates[GLFW_KEY_A]) TranslateRelative(make_float3(-spd, 0, 0));
				if (keystates[GLFW_KEY_D]) TranslateRelative(make_float3(spd, 0, 0));
				if (keystates[GLFW_KEY_W]) TranslateRelative(make_float3(0, 0, spd));
				if (keystates[GLFW_KEY_S]) TranslateRelative(make_float3(0, 0, -spd));
				if (keystates[GLFW_KEY_E]) TranslateRelative(make_float3(0, spd, 0));
				if (keystates[GLFW_KEY_Q]) TranslateRelative(make_float3(0, -spd, 0));

				float2 delta = lastmouse - mousepos;
				float2 rot = dt * 0.00175f * delta;
				TranslateTarget(make_float3(0, rot.y, 0));
				TranslateTarget(make_float3(-rot.x, 0, 0));
			}
			lastmouse = mousepos;
		}

	};

}

namespace vdb
{
	struct ImportingFile
	{
		openvdb::FloatGrid::Ptr grid;
		CoordBBox bbox;

		void Open(const char* file)
		{
			openvdb::io::File newFile(file);
			newFile.open();
			openvdb::GridBase::Ptr baseGrid = newFile.getGrids()->at(0);
			grid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);
			newFile.close();

			bbox = grid->evalActiveVoxelBoundingBox();
		}
	};


	openvdb::FloatGrid::Ptr Transform(openvdb::FloatGrid::Ptr gridIN, float3 trans, float3 rot, float3 scale)
	{
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();

		Vec3d vdbT = Vec3d(trans.x, trans.y, trans.z);
		Vec3d vdbS = Vec3d(scale.x, scale.y, scale.z);
		Vec3d vdbR = Vec3d(rot.x, rot.y, rot.z);

		openvdb::tools::GridTransformer transformer(Vec3s::zero(), vdbS, vdbR, vdbT);
		transformer.transformGrid<openvdb::tools::QuadraticSampler, openvdb::FloatGrid>(*gridIN, *grid);
		grid->tree().prune();
		return grid;
	}

	openvdb::FloatGrid::Ptr Load_cached(const char* file, float3 trans, float3 rot, float scale)
	{
		uint64_t key = (uint64_t)MAPWIDTH << 32 | MAPHEIGHT << 16 | MAPDEPTH;

		std::stringstream stream;
		stream << file << "." << key << ".vdb";
		std::string cachedF(stream.str());
		//std::cout << "trying " << cachedF <<std::endl;

#if 0
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
#endif
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
				grid =  Transform(gridIN, trans, rot, make_float3(scale, scale, scale));

// 				grid = openvdb::FloatGrid::create();
// 
// 				Vec3d trans;
// 				trans.x() = (MAPWIDTH / 2);
// 				trans.z() = (MAPDEPTH / 2);
// 				trans.y() = 0;// (MAPHEIGHT / 2);
// 
// 				openvdb::tools::GridTransformer transformer(Vec3s::zero(), Vec3d(scale, scale, scale), Vec3s::zero(), trans);
// 				transformer.transformGrid<openvdb::tools::QuadraticSampler, openvdb::FloatGrid>(*gridIN, *grid);
// 
// 				grid->tree().prune();

				std::cout << " transform " << t.elapsed() << std::endl;
			}


			openvdb::io::File newFile2(cachedF);
			newFile2.write({ grid });
			newFile2.close();


			return grid;
		}
	}

	void Inject(openvdb::FloatGrid::Ptr grid)
	{
		Timer t;
		struct Injector
		{
			using TreeType = typename openvdb::FloatGrid::TreeType;
			using LeafNode = typename TreeType::ValueOnCIter;
			using IterRange = openvdb::tree::IteratorRange<typename TreeType::ValueOnCIter>;
			void operator()(IterRange& range) const
			{
				for (; range; ++range)
				{
					const openvdb::FloatGrid::ValueOnCIter& iter = range.iterator();
					auto coord = iter.getCoord();
					auto icoord = coord.asVec3I();
					int x = icoord.x();
					int y = icoord.y();
					int z = icoord.z();
					if ((x < MAPWIDTH && x >= 0) && (y < MAPHEIGHT && y >= 0) && (z < MAPDEPTH && z >= 0))
					{
						unsigned char vc = LIGHTBLUE;
						world->Set(x, y, z, vc);
					}
				}
			}
		} proc;
		tbb::parallel_for(Injector::IterRange(grid->tree().cbeginValueOn()), proc);
		std::cout << " inject MT " << t.elapsed() << std::endl;
	}
}


util::camera cam;

void MyGame::Init()
{
	world = GetWorld();
	cam.ct = mat4::LookAt(make_float3(MAPWIDTH, MAPHEIGHT * .66f, -MAPDEPTH * .66f), make_float3(MAPWIDTH / 2, MAPHEIGHT / 2, MAPDEPTH / 2));
	openvdb::initialize();
}
void MyGame::KeyUp(int key)			{cam.keystates[key] = 0;}
void MyGame::KeyDown(int key)		{cam.keystates[key] = 1;}
void MyGame::MouseUp(int button)	{cam.mouseStates[button] = 0;}
void MyGame::MouseDown(int button)	{cam.mouseStates[button] = 1;}
void MyGame::MouseMove(int x, int y){cam.mousepos = make_float2(x, y);}






openvdb::FloatGrid::Ptr gridG;
float3 activeDims;


void MyGame::Tick( float deltaTime )
{

	cam.update(deltaTime);
	world->SetCameraMatrix(cam.ct);
	

	if (ImGui::Begin("OpenVBD importer"))
	{
		ImGui::SetWindowFontScale(2); // dpi scaling for noobs
		if (ImGui::Button("Import..."))
		{
			ImGui::OpenPopup("Import OpenVDB");
		}

		if (ImGui::Button("clear"))
		{
			world->Clear();
		}

		if (ImGui::BeginPopup("Import OpenVDB"))
		{
			ImGui::Text("Import new vdb");

			static char path[MAX_PATH] = "";

			bool changed = ImGui::InputText("File", path, MAX_PATH);
			ImGui::SameLine();
			if (ImGui::Button("..."))
			{
				const char* picked = util::Filepicker();
				if (picked)
				{
					strcpy(path, picked);
					changed = true;
				}
			}

			
			if (FileExists(path))
			{				
				static vdb::ImportingFile file;
				if (changed)
				{
					file.Open(path);
				}

				ImGui::Text("Dimensions: [%d,%d,%d]x[%d,%d,%d] -> [%dx%dx%d]",
					(int)file.bbox.min().x(), (int)file.bbox.min().y(), (int)file.bbox.min().z(),
					(int)file.bbox.max().x(), (int)file.bbox.max().y(), (int)file.bbox.max().z(),
					(int)file.bbox.dim().x(), (int)file.bbox.dim().y(), (int)file.bbox.dim().z());


				CoordBBox newBbox = file.bbox;

				static float3 trans = make_float3(0);
				ImGui::DragFloat3("Translation", (float*)&trans);				

				static float3 rot = make_float3(0);
				ImGui::DragFloat3("Rotation", (float*)&rot);

				static float scale = 1;
				ImGui::SliderFloat("Scale", &scale, 0.01f, 5.0);

				Vec3d bbmin = newBbox.min().asVec3d(), bbmax = newBbox.max().asVec3d();
				bbmin *= (double)scale; bbmax *= (double)scale;
				newBbox.reset(Coord::round(bbmin), Coord::round(bbmax));
				//rotate
				newBbox.translate(Coord(trans.x, trans.y, trans.z));

				
				ImGui::Text("Dimensions: [%d,%d,%d]x[%d,%d,%d] -> [%dx%dx%d]",
					(int)newBbox.min().x(), (int)newBbox.min().y(), (int)newBbox.min().z(),
					(int)newBbox.max().x(), (int)newBbox.max().y(), (int)newBbox.max().z(),
					(int)newBbox.dim().x(), (int)newBbox.dim().y(), (int)newBbox.dim().z());

				if (ImGui::Button("Import"))
				{
					world->Clear();

					auto v = vdb::Load_cached(path, trans, rot, scale);
					vdb::Inject(v);

				}
			}				
			ImGui::EndPopup();
		}
		ImGui::End();
	}

	
// 
// 	static float3 trans = f;
// 	if (trans != f)
// 	{
// 		trans = f;
// 		world->Clear();
// 		gridG = vdb::Transform(gridG, trans, make_float3(0), make_float3(1));
// 		vdb::Inject(gridG, make_float3(0, 0, 0));
// 	}
// 		
	

}