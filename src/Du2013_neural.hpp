#ifndef CELLDU2013_NEURALFROMCELLML_HPP_
#define CELLDU2013_NEURALFROMCELLML_HPP_

//! @file
//!
//! This source file was generated from CellML by chaste_codegen version 0.9.0
//!
//! Model: imtiaz_2002
//!
//! Processed by chaste_codegen: https://github.com/ModellingWebLab/chaste-codegen
//!     (translator: chaste_codegen, model type: Normal)
//! on 2021-09-14 22:23:00
//!
//! <autogenerated>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractStimulusFunction.hpp"
#include "AbstractCardiacCell.hpp"

class CellDu2013_neuralFromCellML : public AbstractCardiacCell
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCardiacCell >(*this);
        
    }

    //
    // Settable parameters and readable variables
    //


public:

    double GetIntracellularCalciumConcentration();
    CellDu2013_neuralFromCellML(boost::shared_ptr<AbstractIvpOdeSolver> pSolver, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);
    ~CellDu2013_neuralFromCellML();
    double GetIIonic(const std::vector<double>* pStateVariables=NULL);
    void EvaluateYDerivatives(double var_chaste_interface__Time__time, const std::vector<double>& rY, std::vector<double>& rDY);

    std::vector<double> ComputeDerivedQuantities(double var_chaste_interface__Time__time, const std::vector<double> & rY);
};

// Needs to be included last
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CellDu2013_neuralFromCellML)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const CellDu2013_neuralFromCellML * t, const unsigned int fileVersion)
        {
            const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
            const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
            ar << p_solver;
            ar << p_stimulus;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, CellDu2013_neuralFromCellML * t, const unsigned int fileVersion)
        {
            boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
            boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
            ar >> p_solver;
            ar >> p_stimulus;
            ::new(t)CellDu2013_neuralFromCellML(p_solver, p_stimulus);
        }

    }

}

#endif // CELLDU2013_NEURALFROMCELLML_HPP_