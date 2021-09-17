#ifndef TESTTRIANGULATEDICCSMC_HPP_
#define TESTTRIANGULATEDICCSMC_HPP_

#include <cxxtest/TestSuite.h>
#include <cassert>

#include <set>

#include "BidomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"
//#include "../src/ICCCBDerivedCa.hpp"
#include "../src/Du2013_neural.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "../src/DummyDerivedCa.hpp"
#include "Debug.hpp"
#include "AbstractElement.hpp"
#include "Node.hpp"
#include "BidomainProblem.hpp"
#include "ChasteEllipsoid.hpp"
#include "AbstractElement.hpp"

using namespace std;

class ICCNwCellFactory : public AbstractCardiacCellFactory<2>
{

private:
    std::set<unsigned> setICCNode;

public:
    ICCNwCellFactory(std::set<unsigned> iccNodes):AbstractCardiacCellFactory<2>(), setICCNode(iccNodes)
    {
        TRACE("Number of node in cell factory: " << setICCNode.size());
    }
    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        unsigned index =  pNode->GetIndex();
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        ChastePoint<2> centre(0.06,0.3);
        ChastePoint<2> radii (0.01,0.01);
        ChasteEllipsoid<2> ellipseRegion(centre, radii);
        ChastePoint<2> myPoint(x, y);
        double eta;
        if( setICCNode.find(index) != setICCNode.end() )
        {
            //ICCCBDerivedCa* cell = new ICCCBDerivedCa(mpSolver, mpZeroStimulus);
            CellDu2013_neuralFromCellML* cell = new CellDu2013_neuralFromCellML(mpSolver, mpZeroStimulus);
            cell->SetParameter("eta", 0.041);
            if(ellipseRegion.DoesContain(myPoint))
            {
                cell->SetParameter("eta", 0.0389);
            }
            return cell;
        }
        return new DummyDerivedCa(mpSolver, mpZeroStimulus);
    }
};

class TestTriangulatedNWFunction : public CxxTest::TestSuite
{
public:
    void TestSimulation() //throw(Exception)
    {
        std::set<unsigned> iccNodes;
        std::vector<ChasteEllipsoid<3> > iccRegion;
        std::vector< c_vector<double,3> > intra_conductivities;
        std::vector< c_vector<double,3> > extra_conductivities;
        DistributedTetrahedralMesh<2,2> mesh;
        unsigned nElements = 0;

        std::string myFile = "MeshNetwork-2D-1601Nodes-3072Elems";
        std::string meshFile = "projects/mesh/ICC2D/" + myFile;
        TrianglesMeshReader<2,2> mesh_reader(meshFile.c_str());

        mesh.ConstructFromMeshReader(mesh_reader);
        nElements = mesh.GetNumLocalElements();
        TRACE("Number of elements: " << nElements);
        double eleIdentify = 0;
        for (DistributedTetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin(); iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            eleIdentify = iter->GetAttribute();
            if (eleIdentify == 1) // ICC=1 and Bath=0
            {
                if(!iter->GetNode(0)->IsBoundaryNode())
                {
                    iccNodes.insert(iter->GetNodeGlobalIndex(0));
                }
                if(!iter->GetNode(1)->IsBoundaryNode())
                {
                    iccNodes.insert(iter->GetNodeGlobalIndex(1));
                }
                if(!iter->GetNode(2)->IsBoundaryNode())
                {
                    iccNodes.insert(iter->GetNodeGlobalIndex(2));
                }
            }
            else
            {

            }
        }

        TRACE("Number of elements: " << nElements);
        TRACE("Number of ICC nodes: " << iccNodes.size());
        TRACE("Total number of nodes: " << mesh.GetNumAllNodes());

        ICCNwCellFactory nwCells(iccNodes);
        BidomainProblem<2> bidomain_problem(&nwCells, true);
        HeartConfig::Instance()->Reset();
	HeartConfig::Instance()->SetSimulationDuration(30000);

        std::string mod = myFile + "-Du_base";
        HeartConfig::Instance()->SetOutputDirectory(mod.c_str());
	HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        std::set<unsigned> tissue_ids;
        tissue_ids.insert(1);
        std::set<unsigned> bath_ids;
        bath_ids.insert(0);
        HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids, bath_ids);
	bidomain_problem.SetMesh( &mesh );
	bidomain_problem.SetWriteInfo();
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.12, 0.12)); // these are quite smaller than cm values
	HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.2, 0.2)); // these are quite smaller than cm values
	HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
	HeartConfig::Instance()->SetCapacitance(2.5);
	HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
	HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, 0.1, 10);
	bidomain_problem.SetWriteInfo();
	bidomain_problem.Initialise();
        HOW_MANY_TIMES_HERE("Check");
	bidomain_problem.Solve();
	HeartEventHandler::Headings();
	HeartEventHandler::Report();
    }
};
#endif
