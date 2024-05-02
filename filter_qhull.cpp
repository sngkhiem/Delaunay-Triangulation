#include "filter_qhull.h"
#include "qhull_tools.h"
#include <vcg/complex/algorithms/convex_hull.h>

using namespace std;
using namespace vcg;

struct MyEdge {
	coordT x1, y1;
	coordT x2, y2;
};

QhullPlugin::QhullPlugin()
{
	typeList = {
		FP_DELAUNAY
	};

	for (ActionIDType tt : types())
		actionList.push_back(new QAction(filterName(tt), this));
}

QhullPlugin::~QhullPlugin()
{
}

QString QhullPlugin::pluginName() const
{
	return "FilterTest2";
}

QString QhullPlugin::filterName(ActionIDType f) const
{
	switch (f) {
	case FP_DELAUNAY: return QString("Delaunay Triangulation 2");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::pythonFilterName(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_DELAUNAY: return QString("Convert cloud point to delaunay triangulation");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::filterInfo(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_DELAUNAY:
		return QString("Convert cloud point to delaunay triangulation");
	default: assert(0);
	}
	return QString("Error: Unknown Filter");
}

QhullPlugin::FilterClass QhullPlugin::getClass(const QAction* a) const
{
	switch (ID(a)) {
	case FP_DELAUNAY: return FilterClass(FilterPlugin::Remeshing);
	default: assert(0);
	}
	return FilterClass(0);
}

RichParameterList QhullPlugin::initParameterList(const QAction* action, const MeshModel& m)
{
	RichParameterList parlst;
	switch (ID(action)) {
	case FP_DELAUNAY:
	default: break; // do not add any parameter for the other filters
	}
	return parlst;
}

coordT *readpointsFromMesh(int* numpoints, int* dimension, MeshModel& m) {
	coordT *points, *coords;

	coords = points = (coordT*) malloc((*numpoints) * (*dimension) * sizeof(coordT));

	int cnt = 0;
	CMeshO::VertexIterator vi;
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
		if (!(*vi).IsD()) {
			for (int ii = 0; ii < *dimension; ++ii)
				*(coords++) = (*vi).P()[ii];
			++cnt;
		}
	assert(cnt == m.cm.vn);

	return (points);
}

bool isInCircumcircle(coordT pX1, coordT pY1, coordT pX2, coordT pY2, coordT pX3, coordT pY3, coordT pX, coordT pY) {
	coordT centerX, centerY;
	coordT a1, b1, d1, a2, b2, d2;
	coordT radius, checkDistance;

	// Make the below calculation easier
	a1 = -2.0 * pX1 + 2.0 * pX2;
	b1 = -2.0 * pY1 + 2.0 * pY2;
	d1 = -pX1 * pX1 - pY1 * pY1 + pX2*pX2 + pY2*pY2;
	a2 = -2.0 * pX1 + 2.0 * pX3;
	b2 = -2.0 * pY1 + 2.0 * pY3;
	d2 = -pX1 * pX1 - pY1 * pY1 + pX3 * pX3 + pY3 * pY3;

	// Using determinant to calculate the coordinate x, and y of center point
	centerX = (d1 * b2 - d2 * b1) / (a1 * b2 - a2 * b1);
	centerY = (d2 * a1 - d1 * a2) / (a1 * b2 - a2 * b1);
	radius  = sqrt((centerX - pX1) * (centerX - pX1) + (centerY - pY1) * (centerY - pY1));

	// Distance between check point and center point
	checkDistance = sqrt((centerX - pX) * (centerX - pX) + (centerY - pY) * (centerY - pY));
	return checkDistance <= radius;
}

// Comapre if two edges are identical
bool edgeCompare(MyEdge edge1, MyEdge edge2) {
	if ((edge1.x1 == edge2.x1 && edge1.y1 == edge2.y1 && edge1.x2 == edge2.x2 &&
		 edge1.y2 == edge2.y2) ||
		(edge1.x1 == edge2.x2 && edge1.y1 == edge2.y2 && edge1.x2 == edge2.x1 &&
		 edge1.y2 == edge2.y1))
		return true;
	return false;
}

std::map<std::string, QVariant> QhullPlugin::applyFilter(
	const QAction*           filter,
	const RichParameterList& par,
	MeshDocument&            md,
	unsigned int& /*postConditionMask*/,
	vcg::CallBackPos* /* cb*/)
{
	qhT  qh_qh = {};
	qhT* qh    = &qh_qh;

	switch (ID(filter)) {
	case FP_DELAUNAY: {
		MeshModel  &m  = *md.mm();
		MeshModel &nm  = *md.addNewMesh("", "Delaunay Triangulation");
		int dim = 3;
		int numpoints = m.cm.vn;
		coordT *points;
		points = readpointsFromMesh(&numpoints, &dim, m);

		/* Create super triangle */
		coordT minX, minY, maxX, maxY;
		minX = minY = 1e9;
		maxX = maxY = -1e9;
		for (int i = 0; i < numpoints; i++)
			for (int j = 0; j < dim; j++) {
				if (j == 0) { // Dimension x
					minX = min(minX, points[i * dim + j]);
					maxX = max(maxX, points[i * dim + j]);
				}
				else if (j == 1) { // Dimension y
					minY = min(minY, points[i * dim + j]);
					maxY = max(maxY, points[i * dim + j]);
				}
			}
		coordT dX = (maxX - minX) * 5.0;
		coordT dY = (maxY - minY) * 5.0;
		// Left-below point
		Point3d p0 = {minX - dX, minY - dY * 3, 0};
		// Left-above point
		Point3d p1 = {minX - dX, maxY + dY * 3, 0};
		// Right-above point
		Point3d p2 = {maxX + dX * 3, maxY + dY * 3, 0};
		// Add super triangle's vertices
		tri::Allocator<CMeshO>::AddVertex(nm.cm, p0);
		tri::Allocator<CMeshO>::AddVertex(nm.cm, p1);
		tri::Allocator<CMeshO>::AddVertex(nm.cm, p2);
		// Add super triangle
		tri::Allocator<CMeshO>::AddFace(nm.cm, p0, p1, p2);
		// Coordinate of super triangle
		coordT superX1 = minX - dX, superY1 = minY - dY * 3;
		coordT superX2 = minX - dX, superY2 = maxY + dY * 3;
		coordT superX3 = maxX + dX * 3, superY3 = maxY + dY * 3;

		/*Triangulation Process*/
		 for (int ii = 0; ii < numpoints; ii++) {
			vector<MyEdge> edges;
			// Add new vertex to current triangulation
			Point3d newVertex = {points[ii * dim], points[ii * dim + 1], 0};
			tri::Allocator<CMeshO>::AddVertex(nm.cm, newVertex);

			// Remove faces with circumcircles containing the vertex
			CMeshO::FaceIterator fi;
			vector<CMeshO::FaceIterator> deleteFaces;
			for (fi = nm.cm.face.begin(); fi != nm.cm.face.end(); fi++) {
				if (!(*fi).IsD()) {
					// Take the coordinate of three vertices of current face
					coordT* curVertices;
					coordT *fpoints, *coords;
					coords = fpoints = (coordT*) malloc(3 * dim * sizeof(coordT));
					for (int i = 0; i < 3; i++)
						for (int j = 0; j < dim; ++j)
							*(coords++) = (*fi).P(i)[j];
					curVertices = fpoints;

					// Check illegal face
					if (isInCircumcircle(
								 curVertices[0],
								 curVertices[1],
								 curVertices[3], curVertices[4],
								 curVertices[6], curVertices[7],
								 points[ii * dim], points[ii * dim + 1])) {
						edges.push_back({curVertices[0], curVertices[1], curVertices[3], curVertices[4]});
						edges.push_back({curVertices[0], curVertices[1], curVertices[6], curVertices[7]});
						edges.push_back({curVertices[3], curVertices[4], curVertices[6], curVertices[7]});
						deleteFaces.push_back(fi);
					}
				}
			}

			// Delete illegal faces
			for (int i = 0; i < deleteFaces.size(); i++)
				tri::Allocator<CMeshO>::DeleteFace(nm.cm, *deleteFaces[i]);

			// Keep the unique edges, because unique edges (which not share with other faces) are boundary edges
			vector<MyEdge> boundaryEdges;
			for (int i = 0; i < edges.size(); i++) {
				bool unique = true;

				for (int j = 0; j < edges.size(); j++)
					if (i != j && edgeCompare(edges[i], edges[j])) {
							unique = false;
							break;
						}

				if (unique)
					boundaryEdges.push_back(edges[i]);
			}

			// Create new faces from the boundary edges
			for (int i = 0; i < boundaryEdges.size(); i++) {
				Point3d p0 = {boundaryEdges[i].x1, boundaryEdges[i].y1, 0};
				Point3d p1 = {boundaryEdges[i].x2, boundaryEdges[i].y2, 0};
				tri::Allocator<CMeshO>::AddFace(nm.cm, p0, p1, newVertex);
			}

			//Clear vector
			edges.clear();
			boundaryEdges.clear();
			deleteFaces.clear();
		}

		 // Remove triangles that share edges with super triangle
		CMeshO::FaceIterator fi;
		vector<CMeshO::FaceIterator> deleteFaces;
		for (fi = nm.cm.face.begin(); fi != nm.cm.face.end(); fi++)
			if (!(*fi).IsD()) {
				// Edges of current face
				coordT* curVertices;
				coordT *fpoints, *coords;
				coords = fpoints = (coordT*) malloc(3 * dim * sizeof(coordT));
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < dim; ++j)
						*(coords++) = (*fi).P(i)[j];
				curVertices = fpoints;

				// Check if 1 of 3 edges of this face is shared with super triangle edges
				// First vertex
				if ((superX1 == curVertices[0] && superY1 == curVertices[1]) ||
					(superX1 == curVertices[3] && superY1 == curVertices[4]) ||
					(superX1 == curVertices[6] && superY1 == curVertices[7]))
					deleteFaces.push_back(fi);
				// Second vertex
				else if ((superX2 == curVertices[0] && superY2 == curVertices[1]) ||
					(superX2 == curVertices[3] && superY2 == curVertices[4]) ||
					(superX2 == curVertices[6] && superY2 == curVertices[7]))
					deleteFaces.push_back(fi);
				// Third vertex
				else if ((superX3 == curVertices[0] && superY3 == curVertices[1]) ||
					(superX3 == curVertices[3] && superY3 == curVertices[4]) ||
					(superX3 == curVertices[6] && superY3 == curVertices[7]))
					deleteFaces.push_back(fi);
				
			}
		for (int i = 0; i < deleteFaces.size(); i++)
			tri::Allocator<CMeshO>::DeleteFace(nm.cm, *deleteFaces[i]);
		log("Succes convert cloud points to triangle mesh (Delaunay Triangulation)");

	} break;
	default: wrongActionCalled(filter);
	}

	return std::map<std::string, QVariant>();
}
MESHLAB_PLUGIN_NAME_EXPORTER(QhullPlugin)
