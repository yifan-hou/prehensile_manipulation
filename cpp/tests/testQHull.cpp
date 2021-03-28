#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/Qhull.h"

#include <iostream>
#include <Eigen/Dense>

using orgQhull::Qhull;
// using orgQhull::QhullError;
// using orgQhull::QhullFacet;
using orgQhull::QhullFacetList;
using orgQhull::QhullPoint;
using orgQhull::RboxPoints;
using orgQhull::QhullVertexList;
using orgQhull::QhullVertexListIterator;
// using orgQhull::QhullVertexSet;

int main() {
    // RboxPoints rbox("100");
    // RboxPoints rbox;
    // rbox.appendPoints("100 D4");
    // Qhull q(rbox, "");
    Eigen::MatrixXd mat(3, 4);
    mat << 0, 2, 0, 0,
           0, 0, 2, 0,
           0, 0, 0, 2;

    double points[] = {
        0, 0, 0,
        2, 0, 0,
        0, 2, 0,
        0, 0, 2
    };
    int dim = 3;
    int num = 4;

    double *data = mat.transpose().data();
    Qhull q("triangle", dim, num, data, "");
    QhullVertexList vertices(q.beginVertex(), q.endVertex());
    std::cout << "Found: " << vertices.count() << " vertices." << std::endl;;

    for (QhullVertexListIterator i = vertices; i.hasNext(); ) {
        std::cout << "id: " << i.peekNext().id() << ", data: ";
        QhullPoint point = i.next().point();
        double *data = point.coordinates();
        for (int j = 0; j < dim; ++j) std::cout << data[j] << ", ";
        std::cout << std::endl;
    }

    std::cout << "Found: " << vertices.count() << " vertices." << std::endl;;

}



