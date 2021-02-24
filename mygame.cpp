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
	
	openvdb::FloatGrid::Ptr Transform(openvdb::FloatGrid::Ptr gridIN, float3 trans, float3 rot, float3 scale)
	{
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();

		Vec3d vdbT = Vec3d(trans.x,trans.y, trans.z);
		Vec3d vdbS = Vec3d(scale.x, scale.y, scale.z);
		Vec3d vdbR = Vec3d(rot.x, rot.y, rot.z);

		openvdb::tools::GridTransformer transformer(Vec3s::zero(), vdbS, vdbR, vdbT);
		transformer.transformGrid<openvdb::tools::QuadraticSampler, openvdb::FloatGrid>(*gridIN, *grid);

		CoordBBox bbox = grid->evalActiveVoxelBoundingBox();
		printf("Dimensions: [%d,%d,%d]x[%d,%d,%d] -> [%dx%dx%d]\n",
			(int)bbox.min().x(), (int)bbox.min().y(), (int)bbox.min().z(),
			(int)bbox.max().x(), (int)bbox.max().y(), (int)bbox.max().z(),
			(int)bbox.dim().x(), (int)bbox.dim().y(), (int)bbox.dim().z());

		grid->tree().prune();
		return grid;
	}

	void Inject(openvdb::FloatGrid::Ptr grid)
	{
		Timer t;
		struct Injector
		{
			openvdb::FloatGrid::Ptr grid;
			float min, max;

			float Remap(float value, float from1, float to1, float from2, float to2) const
			{
				return (value - from1) / (to1 - from1) * (to2 - from2) + from2;
			}

			using TreeType = typename openvdb::FloatGrid::TreeType;
			using LeafNode = typename TreeType::ValueOnCIter;
			using IterRange = openvdb::tree::IteratorRange<typename TreeType::ValueOnCIter>;
			void operator()(IterRange& range) const
			{
				for (; range; ++range)
				{
					const openvdb::FloatGrid::ValueOnCIter& iter = range.iterator();
					auto coord = iter.getCoord();
					Vec3d wcoord = coord.asVec3d();
					Vec3i icoord = Vec3i(wcoord);
					int x = icoord.x();
					int y = icoord.y();
					int z = icoord.z();
					if ((x < MAPWIDTH && x >= 0) && (y < MAPHEIGHT && y >= 0) && (z < MAPDEPTH && z >= 0))
					{
						unsigned char vc = LIGHTBLUE;
						float fval = *iter;
						fval = Remap(fval, min, max, 0, std::numeric_limits<payload>::max());

 						payload p = (payload)fval;

						world->Set(x, y, z, p);
					}
				}
			}
		} proc;

		float min, max;
		grid->evalMinMax(min, max);
		std::cout << " Min val " << min << " Max val "<< max << std::endl;
		RenderParams& p = GetWorld()->GetParams();
		p.dims = make_float2(min, max);

		proc.grid = grid;
		proc.min = min;
		proc.max = max;

		tbb::parallel_for(Injector::IterRange(grid->tree().cbeginValueOn()), proc);
		std::cout << " inject MT " << t.elapsed() << std::endl;
	}

	struct ImportingFile
	{
		
		openvdb::FloatGrid::Ptr grid;
		CoordBBox bbox;	
		float3 trans = make_float3(0);
		float3 rot = make_float3(0);
		float scale = 1;


		void Open(const char* file)
		{
			openvdb::io::File newFile(file);
			newFile.open();
			openvdb::GridBase::Ptr baseGrid = newFile.getGrids()->at(0);
			grid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);
			newFile.close();

			bbox = grid->evalActiveVoxelBoundingBox();
		}

		void Reset()
		{
			grid = nullptr;
			bbox.reset();
		}

		static void ImportUI()
		{
			static ImportingFile* file = nullptr;
			static char path[MAX_PATH] = "";

			if (ImGui::BeginPopup("Import OpenVDB"))
			{
				ImGui::Text("Import new vdb");				

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
					if (changed)
					{
						if (file == nullptr)
						{
							file = new ImportingFile();
						}
						file->Open(path);
					}

					if (file) 
					{
						ImGui::Text("Dimensions: [%d,%d,%d]x[%d,%d,%d] -> [%dx%dx%d]",
							(int)file->bbox.min().x(), (int)file->bbox.min().y(), (int)file->bbox.min().z(),
							(int)file->bbox.max().x(), (int)file->bbox.max().y(), (int)file->bbox.max().z(),
							(int)file->bbox.dim().x(), (int)file->bbox.dim().y(), (int)file->bbox.dim().z());


						CoordBBox newBbox = file->bbox;

						ImGui::DragFloat3("Translation", (float*)&file->trans);
						ImGui::SameLine();
						if (ImGui::Button("Center"))
						{
							Vec3d mi = newBbox.getCenter();
							file->trans = make_float3(MAPWIDTH >> 1, MAPHEIGHT >> 1, MAPDEPTH >> 1) - make_float3((int)mi.x(), (int)mi.y(), (int)mi.z());
						}

						//ImGui::DragFloat3("Rotation", (float*)&rot);
						ImGui::SliderFloat("Scale", &file->scale, 0.01f, 5.f);
						ImGui::SameLine();
						if (ImGui::Button("Fit"))
						{
							float w = newBbox.dim().x();
							float h = newBbox.dim().y();
							float d = newBbox.dim().z();
							w = MAPWIDTH < w ? (float)MAPWIDTH / w : 1;
							h = MAPHEIGHT < h ? (float)MAPHEIGHT / h : 1;
							d = MAPDEPTH < d ? (float)MAPDEPTH / d : 1;
							file->scale = 0.98 * std::min(std::min(w, h), d);
						}

						Vec3d bbmin = newBbox.min().asVec3d(), bbmax = newBbox.max().asVec3d();
						bbmin *= (double)file->scale; bbmax *= (double)file->scale;
						newBbox = CoordBBox(Coord::floor(bbmin), Coord::floor(bbmax));
						//rotate
						newBbox.translate(Coord(file->trans.x, file->trans.y, file->trans.z));


						ImGui::Text("Dimensions: [%d,%d,%d]x[%d,%d,%d] -> [%dx%dx%d]",
							(int)newBbox.min().x(), (int)newBbox.min().y(), (int)newBbox.min().z(),
							(int)newBbox.max().x(), (int)newBbox.max().y(), (int)newBbox.max().z(),
							(int)newBbox.dim().x(), (int)newBbox.dim().y(), (int)newBbox.dim().z());

						if (ImGui::Button("Import"))
						{
							world->Clear();

							auto v = vdb::Transform(file->grid, file->trans, file->rot, make_float3(file->scale));
							vdb::Inject(v);

							delete file;
							file = NULL;
							path[0] = '\0';
						}
					}
				}
				ImGui::EndPopup();
			}
		}
	};

}


util::camera cam;

float3 ld{ 0.3,0.3,0 };

void MyGame::Init()
{
	world = GetWorld();
	cam.ct = mat4::LookAt(make_float3(MAPWIDTH, MAPHEIGHT * .66f, -MAPDEPTH * .66f), make_float3(MAPWIDTH / 2, MAPHEIGHT / 2, MAPDEPTH / 2));
	openvdb::initialize();

	RenderParams& p = GetWorld()->GetParams();
	p.enableSR = 1;
	p.lightDir = normalize(ld);
	p.raycounts = make_int2(16, 4);
	p.gain = 0.2;
	p.scattering = make_float3(1.5, 1.5, 1.5);
	p.absortion = make_float3(0.4, 0.2, 0.1); //1.5
	p.steps = make_float2(1, 3);
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
	printf("%f\n", deltaTime);
	cam.update(deltaTime);
	world->SetCameraMatrix(cam.ct);
	

	if (ImGui::Begin("OpenVBD importer"))
	{
		ImGui::SetWindowFontScale(2); // dpi scaling for noobs
		if (ImGui::Button("Import..."))
		{
			ImGui::OpenPopup("Import OpenVDB");
		}

//		ImGui::ShowDemoWindow();

		RenderParams& p = GetWorld()->GetParams();


		{
			bool sr = p.enableSR == 0 ? false : true;
			ImGui::Checkbox("ShadowRays", &sr);
			p.enableSR = sr ? 1 : 0;
		}

		{
			ImGui::SliderInt2("raycounts", (int*)&p.raycounts,1, 128);
			ImGui::SliderFloat2("steps", (float*)&p.steps, 0, 10);
		}
		ImGui::SliderFloat("gain", &p.gain, 0, 1);

		ImGui::SliderFloat3("scattering", (float*)&p.scattering, 0, 10);
		ImGui::SliderFloat3("absortion", (float*)&p.absortion, 0, 10);


		ImGui::SliderFloat("ld x", &p.lightDir.x, -1, 1);
		ImGui::SliderFloat("ld y", &p.lightDir.y, -1, 1);
		ImGui::SliderFloat("ld z", &p.lightDir.z, -1, 1);
		p.lightDir = normalize(p.lightDir);

		vdb::ImportingFile::ImportUI();

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