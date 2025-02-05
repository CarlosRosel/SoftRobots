/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
*                           Plugin SoftRobots                                 *
*                                                                             *
* This plugin is also distributed under the GNU LGPL (Lesser General          *
* Public License) license with the same conditions than SOFA.                 *
*                                                                             *
* Contributors: Defrost team  (INRIA, University of Lille, CNRS,              *
*               Ecole Centrale de Lille)                                      *
*                                                                             *
* Contact information: https://project.inria.fr/softrobot/contact/            *
******************************************************************************/
#pragma once

#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/logging/Messaging.h>

#include <SoftRobots/component/constraint/model/MagnetModel.h>
#include <Eigen/Dense>
//#include <Eigen/Geometry>  // Para cuaterniones
//using namespace std;
using namespace Eigen;

namespace softrobots::constraint
{

using sofa::core::objectmodel::ComponentState;
using sofa::core::VecCoordId;
using sofa::core::ConstVecCoordId ;
using sofa::helper::WriteAccessor ;
using sofa::helper::ReadAccessor ;
using sofa::type::vector ;
using sofa::type::Vec;
using sofa::type::Vec3;
using sofa::type::RGBAColor;



template<class DataTypes>
MagnetModel<DataTypes>::MagnetModel(MechanicalState* object)
    : Inherit1(object)
    , d_indices(initData(&d_indices, "indices",
                                 "If indices size is lower than target size, \n"
                                 "some target will not be considered"))

    , d_weight(initData(&d_weight, sofa::type::vector<Real>(Deriv::total_size, 1.), "weight",
                          "The parameter sets a weight to the minimization."))

    , d_directions(initData(&d_directions,"directions",
                          "The parameter directions allows to specify the directions in \n"
                          "which you want to solve the position."))

    , d_Jacobian(initData(&d_Jacobian,"JacobianBq",
                          "Jacobian relating node motion to B change\n"
                          "."))

    , d_useDirections(initData(&d_useDirections,"useDirections",
                              "The parameter useDirections allows to select the directions in \n"
                              "which you want to solve the position. If unspecified, the default \n"
                              "values are all true."))

    , d_delta(initData(&d_delta, "delta","Distance to target"))
{
    d_delta.setReadOnly(true);

    this->addUpdateCallback("updateWeight", {&d_weight}, [this](const sofa::core::DataTracker& t)
                            {
                                SOFA_UNUSED(t);
                                auto weight = sofa::helper::getWriteAccessor(d_weight);
                                if (weight.size() != Deriv::total_size)
                                {
                                    msg_info() << "Wrong size for the data field weight, " << weight.size() <<
                                        " instead of " << Deriv::total_size << ". Resizing, with weight[0] as the default value.";
                                    Real w = weight.empty()? 1.: weight[0];
                                    d_weight.setValue(sofa::type::vector<Real>(Deriv::total_size, w));
                                }
                                return sofa::core::objectmodel::ComponentState::Valid;
                            }, {});
}


template<class DataTypes>
MagnetModel<DataTypes>::~MagnetModel()
{
}


template<class DataTypes>
void MagnetModel<DataTypes>::init()
{
    d_componentState = ComponentState::Valid;
    Inherit1::init();

    if(m_state==nullptr)
    {
        msg_error() << "There is no mechanical state associated with this node. "
                        "the object is deactivated. "
                        "To remove this error message fix your scene possibly by "
                        "adding a MechanicalObject." ;
        d_componentState = ComponentState::Invalid;
        return;
    }

    internalInit();
}


template<class DataTypes>
void MagnetModel<DataTypes>::reinit()
{
    internalInit();
}


template<class DataTypes>
void MagnetModel<DataTypes>::internalInit()
{
    if(!d_directions.isSet())
        setDefaultDirections();
    else
        normalizeDirections();

    if(!d_useDirections.isSet())
    {
        setDefaultUseDirections();
    }
    else
    {
        const auto useDirections = sofa::helper::getReadAccessor(d_useDirections);
        if (std::find(useDirections.begin(), useDirections.end(), true) == useDirections.end())
        {
            setDefaultUseDirections();
            msg_warning(this) << "No direction given in useDirection. Set default all.";
        }
    }

    if(!d_indices.isSet())
    {
        msg_warning(this) <<"Indices not defined. Default value assigned 0.";
        setIndicesDefaultValue();
    }

    if(d_indices.getValue().size() > m_state->getSize())
    {
        msg_warning(this) <<"Indices size can not be larger than the number of point in the context. Launch resize process.";
        resizeIndicesRegardingState();
    }

    if(d_indices.getValue().size() == 0)
    {
        msg_error(this) <<"Indices size is zero. The component will not work.";
        d_componentState = ComponentState::Invalid;
        return;
    }

    checkIndicesRegardingState();
}


template<class DataTypes>
void MagnetModel<DataTypes>::checkIndicesRegardingState()
{
    ReadAccessor<sofa::Data<VecCoord> > positions = m_state->readPositions();

    if(d_indices.getValue().size() > positions.size())
    {
        msg_error(this) << "Indices size is larger than mechanicalState size" ;
        d_componentState = ComponentState::Invalid;
        return;
    }

    const auto& indices = d_indices.getValue();
    for(unsigned int i=0; i<indices .size(); i++)
    {
        if (positions.size() <= indices[i])
        {
            msg_error(this) << "Index at index " << i << " is too large regarding mechanicalState [position] size" ;
            d_componentState = ComponentState::Invalid;
            return;
        }
    }
}




template<class DataTypes>
void MagnetModel<DataTypes>::setIndicesDefaultValue()
{
    WriteAccessor<sofa::Data<vector<unsigned int> > > defaultIndices = d_indices;
    defaultIndices.resize(1);
}

template<class DataTypes>
void MagnetModel<DataTypes>::resizeIndicesRegardingState()
{
    WriteAccessor<sofa::Data<vector<unsigned int>>> indices = d_indices;
    indices.resize(m_state->getSize());
}


Eigen::Vector3d Calculo_B_Test(double x,double y,double z,double mu_mag_delGrafico, double mu_hat_x,double mu_hat_y,double mu_hat_z)
{
    // Definición de variables
    Vector3d Distancia_r(x, y, z);
    //std::cerr << "Valor de Distancia_r: " << Distancia_r << std::endl;
    double length_r = Distancia_r.norm();  // Magnitud del vector
    //std::cerr << "Valor de length_r: " << length_r << std::endl;

    Vector3d r_hat = Distancia_r / length_r; // Vector unitario
    // cout << "r_hat: " << r_hat.transpose() << endl;

    Vector3d Mu_hat(mu_hat_x, mu_hat_y,mu_hat_z);  // Dirección del momento magnético (vector unitario)
    // cout << "Mu_hat : " << Mu_hat.transpose() << endl;



    // Multiplicación por mu (magnitud)
    Eigen::Vector3d mu = Mu_hat * mu_mag_delGrafico;

    // Producto tensorial r_hat * r_hat^T
    Eigen::Matrix3d AAA = 3 * (r_hat * r_hat.transpose()) - Eigen::Matrix3d::Identity();

    // Multiplicamos AAA por Mu_hat y luego por mu
    Eigen::Vector3d numerador = AAA * mu;
    double denominador = 4 * M_PI * std::pow(std::abs(length_r), 3);

    // Campo magnético B (vector)
    Eigen::Vector3d Campo_Magnetico_resultado = numerador / denominador;
    Campo_Magnetico_resultado *= 1e15;  // Multiplicamos por 10^12 (en unidades apropiadas)

    // cout << "Campo_Magnetico_resultado: " << Campo_Magnetico_resultado.transpose() << endl;


    return Campo_Magnetico_resultado;
}
template<class DataTypes>
void MagnetModel<DataTypes>::buildConstraintMatrix(const ConstraintParams* cParams,
                                                        DataMatrixDeriv &cMatrix,
                                                        unsigned int &cIndex,
                                                        const DataVecCoord &x)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return;

    SOFA_UNUSED(cParams);
    SOFA_UNUSED(x);



// Aqui implementar derivadas


    d_constraintIndex.setValue(cIndex);
    const auto& constraintIndex = sofa::helper::getReadAccessor(d_constraintIndex);
    MatrixDeriv& column = *cMatrix.beginEdit();
    const auto& indices = sofa::helper::getReadAccessor(d_indices);
    sofa::Index sizeIndices = indices.size();
    const auto& useDirections = sofa::helper::getReadAccessor(d_useDirections);
    auto directions = sofa::helper::getWriteAccessor(d_directions);
    auto Jacobian = sofa::helper::getWriteAccessor(d_Jacobian);
    auto weight = sofa::helper::getReadAccessor(d_weight);

//    std::cout << "datos d epos AAAAAAA: " << x << std::endl;
    auto& data = x;  // Usamos directamente el ReadAccessor    double dBx_drz = 0;
    Eigen::Vector3d B_calculada;  // Variable global o de ámbito extendido
    auto readAccessor = sofa::helper::getReadAccessor(data);  // Obtén un ReadAccessor
    Eigen::Vector3d B_0;  // Variable global o de ámbito extendido
    Eigen::Vector3d B_1;  // Variable global o de ámbito extendido
    Eigen::Vector3d B_2;  // Variable global o de ámbito extendido
    Eigen::Vector3d B_3;  // Variable global o de ámbito extendido
    Eigen::Vector3d dBdR_z;  // Variable global o de ámbito extendido
    Eigen::Vector3d dBdR_x;  // Variable global o de ámbito extendido
    Eigen::Vector3d dBdR_y;  // Variable global o de ámbito extendido
    const auto Cambio = 0.01;



    // Ahora puedes iterar sobre los datos subyacentes
    for (const auto& coord : readAccessor) {
            // Acceder a todas las componentes de cada Vec<2, double> y mostrar las tres componentes si es posible
            // msg_warning() << "Coord: (" << coord[0] << ", " << coord[1] << ", " << coord[2] << ")";  // Asumiendo que 'Vec<2, double>' tiene tres componentes
            Eigen::Quaterniond MiR(coord[6], coord[3], coord[4], coord[5]);  // (w, x, y, z)
            Eigen::Matrix3d rotation_matrix = MiR.toRotationMatrix();
            // std::cout << "Matriz de rotación:\n" << rotation_matrix << std::endl;
            // Extraer la última columna
            const double mu_x = rotation_matrix(0, 2);  // Elemento (0, 2)
            const double mu_y = rotation_matrix(1, 2);  // Elemento (1, 2)
            const double mu_z = rotation_matrix(2, 2);  // Elemento (2, 2)

            const double ajuste_x =3.75 ;
            const double ajuste_y = -5;
            const double ajuste_z = 2.7;
    //        const double ajuste_z = 2.7 - 15; // Este ajuste era para corroborar calculo en c++ con el de python
            B_calculada = Calculo_B_Test(coord[0]- ajuste_x,coord[1]- ajuste_y,coord[2],3.81e-9,mu_x,mu_y,mu_z);
//            std::cout << "Campo magnético B_c alculado c++ AAAAAAAAAAAAAA: " << B_calculada.transpose() << std::endl;
            B_0 = Calculo_B_Test(coord[0]- ajuste_x,coord[1]- ajuste_y,coord[2],3.81e-9,mu_x,mu_y,mu_z);

            // dB/drX
            B_2 = Calculo_B_Test(coord[0] - ajuste_x + Cambio,coord[1]- ajuste_y,coord[2],3.81e-9,mu_x,mu_y,mu_z);
            dBdR_x = (B_2 - B_0)/Cambio;
            std::cout << "dBdR_x" << dBdR_x.transpose() << std::endl;
            // dB/drY
            B_3 = Calculo_B_Test(coord[0] - ajuste_x,coord[1]- ajuste_y+ Cambio,coord[2],3.81e-9,mu_x,mu_y,mu_z);
            dBdR_y = (B_3 - B_0)/Cambio;
            std::cout << "dBdR_y" << dBdR_y.transpose() << std::endl;
            // dB/drZ
            B_1 = Calculo_B_Test(coord[0] - ajuste_x,coord[1]- ajuste_y,coord[2]    + Cambio    ,3.81e-9,mu_x,mu_y,mu_z);
            dBdR_z = (B_1 - B_0)/Cambio;
            std::cout << "dBdR_z" << dBdR_z.transpose() << std::endl;

        }


    double dBx_drx = dBdR_x[0];
    double dBy_drx = dBdR_x[1];
    double dBz_drx = dBdR_x[2];
    double dBx_dry = dBdR_y[0];
    double dBy_dry = dBdR_y[1];
    double dBz_dry = dBdR_y[2];
    double dBx_drz = dBdR_z[0];
    double dBy_drz = dBdR_z[1];
    double dBz_drz = dBdR_z[2];



//    Jacobian.push_back(VecDeriv(0.0, 0.0, dBx_drz, 0.0, 0.0, 0.0));


//    Jacobian[0] = sofa::type::Vec<6, double>(17114, 0.0, 0.0, 0.0, 0.0, 0.0);  // Asignar un nuevo valor
//    Jacobian[1] = sofa::type::Vec<6, double>(0.0, 17114, 0.0, 0.0, 0.0, 0.0);  // Asignar un nuevo valor
//    Jacobian[2] = sofa::type::Vec<6, double>(-126.7, -126.7, -33978, 0.0, 0.0, 0.0);  // Asignar un nuevo valor

        Jacobian[0] = sofa::type::Vec<6, double>(dBx_drx, dBx_dry, dBx_drz, 0.0, 0.0, 0.0);  // Asignar un nuevo valor
        Jacobian[1] = sofa::type::Vec<6, double>(dBy_drx, dBy_dry, dBy_drz, 0.0, 0.0, 0.0);  // Asignar un nuevo valor
        Jacobian[2] = sofa::type::Vec<6, double>(dBz_drx, dBz_dry, dBz_drz, 0.0, 0.0, 0.0);  // Asignar un nuevo valor

    unsigned int index = 0;

    for (unsigned j = 0; j < 3; j++) {
        MatrixDerivRowIterator rowIterator = column.writeLine(constraintIndex + index);
//        std::cout << "datos de directions: " << directions[j] << std::endl;

        rowIterator.setCol(0, Jacobian[j]);
        index++;
    }


    // for (unsigned i=0; i<sizeIndices; i++)
    // {
    //     for(sofa::Size j=0; j<Deriv::total_size; j++)
    //     {
    //         if(useDirections[j])
    //         {
    //             MatrixDerivRowIterator rowIterator = column.writeLine(constraintIndex+index);
    //             std::cout << "Tipo de directions: " << typeid(directions).name() << std::endl;
    //             rowIterator.setCol(indices[i], directions[j]);
    //             index++;
    //         }
    //     }
    // }

    cIndex += index;
    cMatrix.endEdit();
    
    m_nbLines = cIndex - constraintIndex;
}


template<class DataTypes>


void MagnetModel<DataTypes>::storeResults(vector<double> &delta)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return;

    d_delta.setValue(delta);
}


template<class DataTypes>
void MagnetModel<DataTypes>::setDefaultDirections()
{
    VecDeriv directions(Deriv::total_size);
    for(sofa::Size i=0; i<Deriv::total_size; i++)
        directions[i][i] = 1.;
    d_directions.setValue(directions);
    d_Jacobian.setValue(directions);
}

template<class DataTypes>
void MagnetModel<DataTypes>::setDefaultUseDirections()
{
    Vec<Deriv::total_size, bool> useDirections;
    useDirections.assign(true);
    d_useDirections.setValue(useDirections);
}


template<class DataTypes>
void MagnetModel<DataTypes>::normalizeDirections()
{

    WriteAccessor<sofa::Data<VecDeriv>> directions = d_directions;
    directions.resize(Deriv::total_size);
    for(unsigned int i=0; i<Deriv::total_size; i++)
        directions[i].normalize();
}


template<class DataTypes>
void MagnetModel<DataTypes>::draw(const VisualParams* vparams)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return;

    if (!vparams->displayFlags().getShowInteractionForceFields())
        return;

    ReadAccessor<sofa::Data<VecCoord> > positions = m_state->readPositions();
    ReadAccessor<sofa::Data<sofa::type::vector<sofa::Index>> > indices = d_indices;
    vector<Coord> points;
    points.reserve(indices.size());
    for (unsigned int i=0; i<indices.size(); i++)
    {
        points.push_back(positions[indices[i]]);
    }
    drawPoints(vparams, points, 10.0f, RGBAColor::green());
}

} // namespace

