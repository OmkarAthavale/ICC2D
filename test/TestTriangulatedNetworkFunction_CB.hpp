#ifndef TESTTRIANGULATEDICCSMC_HPP_
#define TESTTRIANGULATEDICCSMC_HPP_

#include <cxxtest/TestSuite.h>
#include <cassert>

#include <set>

#include "BidomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "../src/ICCCBDerivedCa.hpp"
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
        ChastePoint<2> centre(0.02,0.02);
        ChastePoint<2> radii (0.01, 0.01);
        ChasteEllipsoid<2> ellipseRegion(centre, radii);
        ChastePoint<2> myPoint(x, y);

        if( setICCNode.find(index) != setICCNode.end() )
        {
            ICCCBDerivedCa* cell = new ICCCBDerivedCa(mpSolver, mpZeroStimulus);
            cell->SetParameter("ode_time_step",0.1);
            cell->SetParameter("IP3Par", 0.0006);
            cell->SetParameter("V_excitation", -65);

            //if(ellipseRegion.DoesContain(myPoint))
            //    cell->SetParameter("V_excitation", -68);

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

        //int argc = *(CommandLineArguments::Instance()->p_argc);
		//char** argv = *(CommandLineArguments::Instance()->p_argv);

		//std::string myFile(argv[1]);

        std::set<unsigned> iccNodes;
        std::vector<ChasteEllipsoid<3> > iccRegion;
        std::vector< c_vector<double,3> > intra_conductivities;
        std::vector< c_vector<double,3> > extra_conductivities;
        DistributedTetrahedralMesh<2,2> mesh;
        unsigned nElements = 0;

        //std::string myFile = "63552ICCGHuRx20z08stacktile3x31_2_B";
        //std::string meshFile = "projects/mesh/ICC2D/" + myFile + ".1";
        std::string myFile = "MeshNetwork-2D-1601Nodes-3072Elems";
        std::string meshFile = "projects/mesh/ICC2D/" + myFile;

        //std::string meshFile = "projects/commandLineV2/test/network/" + myFile + ".3";
        TrianglesMeshReader<2,2> mesh_reader(meshFile.c_str());

        mesh.ConstructFromMeshReader(mesh_reader);
        nElements = mesh.GetNumLocalElements();
        TRACE("Number of elements: " << nElements);
        double eleIdentify = 0;
        c_matrix<double, 2, 2> matJac;
        c_matrix<double, 2, 2> matInvJac;
        c_vector<double, 3> circum;
        double determinant = 0;
        for (DistributedTetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin(); iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            //Element<2,2>* presentEle = *iter;
            eleIdentify = iter->GetAttribute();
            //Node<2>* nodeInfo = 0;
            if (eleIdentify > 0.2)
            {
                iter->CalculateJacobian(matJac, determinant);
                iter->CalculateInverseJacobian(matJac, determinant, matInvJac);
                circum = iter->CalculateCircumsphere(matJac, matInvJac);
                ChastePoint<3> centre(circum(0),circum(1), 0);
                double radius = sqrt(circum(2));
                ChastePoint<3> radii (1.1 * radius, 1.1 * radius, 0);

                ChasteEllipsoid<3> ellipseRegion(centre, radii);
                iccRegion.push_back(ellipseRegion);
                intra_conductivities.push_back( Create_c_vector(0.12, 0.12, 0.0));
                extra_conductivities.push_back( Create_c_vector(0.02, 0.02, 0.0));
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
        TRACE("Number of hetero area: " << iccRegion.size());

        ICCNwCellFactory nwCells(iccNodes);
        BidomainProblem<2> bidomain_problem(&nwCells, true);
        HeartConfig::Instance()->Reset();
	      HeartConfig::Instance()->SetSimulationDuration(1);

        std::string mod = myFile + "-CB-v0";
        HeartConfig::Instance()->SetOutputDirectory(mod.c_str());
	      HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        std::set<unsigned> tissue_ids;
        //tissue_ids.insert(0);
        tissue_ids.insert(1);
        std::set<unsigned> bath_ids;
        bath_ids.insert(0);
        //bath_ids.insert(3);
        //bath_ids.insert(4);
        HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids, bath_ids);
	      bidomain_problem.SetMesh( &mesh );
	      bidomain_problem.SetWriteInfo();
        HeartConfig::Instance()->SetConductivityHeterogeneitiesEllipsoid(iccRegion, intra_conductivities, extra_conductivities);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.012, 0.012));
	      HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.002, 0.002));
	      HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
	      HeartConfig::Instance()->SetCapacitance(2.5);
	      HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
	      HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, 0.1, 1);
	      bidomain_problem.SetWriteInfo();
	      bidomain_problem.Initialise();
        HOW_MANY_TIMES_HERE("Check");
	      bidomain_problem.Solve();
	      HeartEventHandler::Headings();
	      HeartEventHandler::Report();
    }
};

#endif
