

#ifndef _MESHFRAME_STATIC_TET_MESH_H_
#define _MESHFRAME_STATIC_TET_MESH_H_

#include "BaseTMesh.h"
#include "TetMeshTypeDefs.h"
#include "../Types/TypeDefs.h"
#include "HoudiniGeoIO.h"
#include <vector>

namespace MF
{
	namespace TetMesh
	{

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		class CTMeshStatic : public CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>
		{
		public:
			typedef CTMeshStatic<DType, TVType, VType, HEType, TEType, EType, HFType, FType, TType>* Ptr;
			typedef std::shared_ptr<CTMeshStatic<DType, TVType, VType, HEType, TEType, EType, HFType, FType, TType>> SharedPtr;

			/*!
				return all the vertex position as a 3xN matrix
			*/
			TVerticesMat<DType>& vertPos() { return mVertPos; };
			TTetIdsMat& tetVIds() { return mTetVIds; };

			TVec3Block<DType> vert(size_t vId) { return mVertPos.block<3, 1>(0, vId); }
			Vec4BlockI tet(size_t tId) { return mTetVIds.block<4, 1>(0, tId);  }

			/*! access the vertex with ID */
			virtual VertexType* idVertex(int id) { return &mVContainer[id]; };

			/*! access the tet with ID */
			virtual TetType* idTet(int id) { return &mTContainer[id]; };

			/*!
			Load tet mesh from a ".t" file
			*/
			void load_t(const char* input, bool checkOrientation = false);

			/*!
			Load tet mesh from a ".geo" file
			*/
			void load_geo(const char* input, bool checkOrientation = false);

		protected:
			TVerticesMat<DType> mVertPos;
			TTetIdsMat mTetVIds;

		};

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		void CTMeshStatic<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::load_t(const char* input, bool checkOrientation)
		{
			std::ostringstream oss;
			std::istringstream iss;

			addVProp(mVHFArrayHandle);
			addVProp(mVTEArrayHandle);

			m_maxVertexId = 0;

			std::fstream is(input, std::fstream::in);

			if (is.fail())
			{
				fprintf(stderr, "Error in opening file %s\n", input);
				return;
			}

			char buffer[MAX_LINE];

			m_nVertices = 0;
			m_nTets = 0;
			m_nEdges = 0;

			while (!is.eof())
			{
				is.getline(buffer, MAX_LINE);
				std::string line(buffer);
				line = strutilTetMesh::trim(line);
				strutilTetMesh::Tokenizer stokenizer(line, " \r\n");

				stokenizer.nextToken();
				std::string token = stokenizer.getToken();

				if (token == "Vertex") m_nVertices++;
				if (token == "Tet") m_nTets++;
				if (token == "Edge") m_nEdges++;
			}

			is.clear();              // forget we hit the end of file
			is.seekg(0, std::ios::beg);   // move to the start of the file

			mVertPos.resize(3, m_nVertices);
			mTetVIds.resize(4, m_nTets);

			//read in the vertices
			int vId = 0;
			std::map<int, int> tFileIdToVId;
			for (int i = 0; i < m_nVertices && is.getline(buffer, MAX_LINE); i++)
			{
				std::string line(buffer);
				line = strutilTetMesh::trim(line);
				strutilTetMesh::Tokenizer stokenizer(line, " \r\n");

				stokenizer.nextToken();
				std::string token = stokenizer.getToken();

				if (token != "Vertex")
				{
					fprintf(stderr, "File Format Error\r\n");
					return;
				}

				stokenizer.nextToken();
				token = stokenizer.getToken();
				int vIdTFile = strutilTetMesh::parseString<int>(token, iss);

				TVec3<DType> p;
				for (int k = 0; k < 3; k++)
				{
					stokenizer.nextToken();
					std::string token = stokenizer.getToken();
					p[k] = strutilTetMesh::parseString<float>(token, iss);
				}

				VertexType* v = createVertexWithIndex();
				v->id() = vId;
				v->setPVertPos(&mVertPos);
				v->position() = p;

				tFileIdToVId.insert({ vIdTFile, vId });
				++vId;
				if (!stokenizer.nextToken("\t\r\n")) continue;
				token = stokenizer.getToken();
			}

			//read in tets 
			int tid = 0;
			for (int id = 0; id < m_nTets && is.getline(buffer, MAX_LINE); id++)
			{
				int vIds[4];

				std::string line(buffer);
				line = strutilTetMesh::trim(line);
				strutilTetMesh::Tokenizer stokenizer(line, " \r\n");

				stokenizer.nextToken();
				std::string token = stokenizer.getToken();

				if (token != "Tet")
				{
					fprintf(stderr, "File Format Error\r\n");
					return;
				}

				//skip the first "4" in the line
				stokenizer.nextToken();
				token = stokenizer.getToken();
				// int tid = strutilTetMesh::parseString<int>(token, iss);


				for (int k = 0; k < 4; k++)
				{
					stokenizer.nextToken();
					std::string token = stokenizer.getToken();
					int vIdTFile = strutilTetMesh::parseString<int>(token, iss);
					vIds[k] = tFileIdToVId[vIdTFile];
				}

				TetType* pT = createTetWithIndex();
				pT->id() = tid;

				if (checkOrientation) {
					_construct_tet_orientation(pT, tid, vIds);
				}
				else {
					_construct_tet(pT, tid, vIds);
				}
				mTetVIds.block<4, 1>(0, tid) << vIds[0], vIds[1], vIds[2], vIds[3];
				tid++;

				if (!stokenizer.nextToken("\t\r\n")) continue;
		
			}

			_construct_faces();
			_construct_edges();

			for (int id = 0; id < m_nEdges && is.getline(buffer, MAX_LINE); id++)
			{
				std::string line(buffer);
				line = strutilTetMesh::trim(line);
				strutilTetMesh::Tokenizer stokenizer(line, " \r\n");

				stokenizer.nextToken();
				std::string token = stokenizer.getToken();

				if (token != "Edge")
				{
					fprintf(stderr, "File Format Error\r\n");
					return;
				}

				stokenizer.nextToken();
				token = stokenizer.getToken();
				int id1 = strutilTetMesh::parseString<int>(token, iss);

				stokenizer.nextToken();
				token = stokenizer.getToken();
				int id2 = strutilTetMesh::parseString<int>(token, iss);

				VertexType* pV1 = idVertex(id1);
				VertexType* pV2 = idVertex(id2);

				EdgeType* pE = VertexEdge(pV1, pV2);

				if (!stokenizer.nextToken("\t\r\n"))
				{
					continue;
				}

				token = stokenizer.getToken();
			}

			m_nEdges = (int)mEContainer.size();

			is.close();

			for (auto vIter = mVContainer.begin(); vIter != mVContainer.end(); vIter++)
			{
				VertexType* pV = *vIter;
				if (pV->id() > m_maxVertexId)
				{
					m_maxVertexId = pV->id();
				}
			}

			// label the boundary for faces and vertices
			for (auto fIter = mFContainer.begin(); fIter != mFContainer.end(); ++fIter)
			{
				FaceType* pF = *fIter;
				if (this->FaceLeftHalfFace(pF) == NULL || this->FaceRightHalfFace(pF) == NULL)
				{
					pF->boundary() = true;
					HalfFaceType* pH =
						FaceLeftHalfFace(pF) == NULL ? FaceRightHalfFace(pF) : FaceLeftHalfFace(pF);
					//added by Anka, mark edge as boundary
					HalfEdgeType* pHE = (HalfEdgeType*)pH->half_edge();

					for (int i = 0; i < 3; ++i)
					{
						EdgeType* pE = HalfEdgeEdge(pHE);
						int vid = pH->key(i);
						VertexType* v = idVertex(vid);
						v->boundary() = true;
						pE->boundary() = true;
						pHE = HalfEdgeNext(pHE);
					}
				}
			}

			for (auto vIter = mVContainer.begin(); vIter != mVContainer.end(); vIter++)
			{
				VertexType* pV = *vIter;
				pV->edges()->shrink_to_fit();
				pV->tvertices()->shrink_to_fit();
			}

			removeVProp(mVTEArrayHandle);

		}


		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		void CTMeshStatic<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::load_geo(const char* input, bool checkOrientation)
		{
			// 已有代码
			std::string p = input;
			HoudiniGeoIO geo(p);
			auto indices = geo.getIndices();
			auto positions = geo.getPositions();
			// 1. 添加临时属性
			addVProp(mVHFArrayHandle);
			addVProp(mVTEArrayHandle);

			// 2. 初始化网格参数
			m_maxVertexId = 0;
			m_nVertices = (int)(positions.size() / 3);  // positions是一维数组，每3个float表示一个顶点
			m_nTets = (int)(indices.size() / 4);        // indices是一维数组，每4个int表示一个四面体
			m_nEdges = 0;                        // 边数后面计算

			// 3. 调整矩阵大小
			mVertPos.resize(3, m_nVertices);
			mTetVIds.resize(4, m_nTets);

			// 4. 创建顶点
			for (int vId = 0; vId < m_nVertices; vId++) {
				TVec3<DType> p;
				p[0] = positions[vId * 3];
				p[1] = positions[vId * 3 + 1];
				p[2] = positions[vId * 3 + 2];

				VertexType* v = createVertexWithIndex();
				v->id() = vId;
				v->setPVertPos(&mVertPos);
				v->position() = p;
			}

			// 5. 创建四面体
			for (int tid = 0; tid < m_nTets; tid++) {
				int vIds[4];
				for (int k = 0; k < 4; k++) {
					vIds[k] = indices[tid * 4 + k];
				}

				TetType* pT = createTetWithIndex();
				pT->id() = tid;

				if (checkOrientation) {
					_construct_tet_orientation(pT, tid, vIds);
				}
				else {
					_construct_tet(pT, tid, vIds);
				}

				mTetVIds.block<4, 1>(0, tid) << vIds[0], vIds[1], vIds[2], vIds[3];
			}
			

			// 6. 构建面和边
			_construct_faces();
			_construct_edges();
			m_nEdges = (int)mEContainer.size();
			
			// 7. 设置最大顶点ID
			for (auto vIter = mVContainer.begin(); vIter != mVContainer.end(); vIter++) {
				VertexType* pV = *vIter;
				if (pV->id() > m_maxVertexId) {
					m_maxVertexId = pV->id();
				}
			}

			// 8. 标记边界面和顶点
			for (auto fIter = mFContainer.begin(); fIter != mFContainer.end(); ++fIter) {
				FaceType* pF = *fIter;
				if (this->FaceLeftHalfFace(pF) == NULL || this->FaceRightHalfFace(pF) == NULL) {
					pF->boundary() = true;
					HalfFaceType* pH = 
						FaceLeftHalfFace(pF) == NULL ? FaceRightHalfFace(pF) : FaceLeftHalfFace(pF);
					HalfEdgeType* pHE = (HalfEdgeType*)pH->half_edge();

					for (int i = 0; i < 3; ++i) {
						EdgeType* pE = HalfEdgeEdge(pHE);
						int vid = pH->key(i);
						VertexType* v = idVertex(vid);
						v->boundary() = true;
						pE->boundary() = true;
						pHE = HalfEdgeNext(pHE);
					}
				}
			}

			// 9. 优化内存
			for (auto vIter = mVContainer.begin(); vIter != mVContainer.end(); vIter++) {
				VertexType* pV = *vIter;
				pV->edges()->shrink_to_fit();
				pV->tvertices()->shrink_to_fit();
			}

			// 10. 移除临时属性
			removeVProp(mVTEArrayHandle);
		}
	}
}

#endif