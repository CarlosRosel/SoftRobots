
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

    , d_PosSensor(initData(&d_PosSensor,"PosSensor",
                          "PosSensor Posicion del sensor\n"
                          "."))
    , d_mum(initData(&d_mum,"mum",
                          "mum \n"
                          "."))
    , d_JacobianFactor(initData(&d_JacobianFactor,"JacobianFactor",
                          "JacobianFactor \n"
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
    if(!d_PosSensor.isSet())
    {
        // setDefaultDirections();
    }
    else
    {
        const auto PosSensor = sofa::helper::getReadAccessor(d_PosSensor);
        // std::cout << "ULTIMO POSSENSOR---------------------------: " << PosSensor << std::endl;
        // std::cout << "ULTIMO POSSENSOR---------------------------0: " << PosSensor[0][0] << std::endl;
    }
    if(!d_mum.isSet())
    {
        // setDefaultDirections();
    }

    else
    {
        const auto mum = sofa::helper::getReadAccessor(d_mum);
        std::cout << "ULTIMO mum---------------------------: " << mum[0][0] << std::endl;
    }

    if(!d_JacobianFactor.isSet() || d_JacobianFactor.getValue().empty())
        {
            typename DataTypes::VecDeriv defaultFactor(1);
            for(size_t i = 0; i < DataTypes::Deriv::total_size; ++i)
                defaultFactor[0][i] = 1.0;
            
            d_JacobianFactor.setValue(defaultFactor);
            std::cout << "ULTIMO JacobianFactor---------------------------: " << d_JacobianFactor.getValue()[0][0] << std::endl;
        }

    else
    {
         const auto JacobianFactor = sofa::helper::getReadAccessor(d_JacobianFactor);
        std::cout << "ULTIMO JacobianFactor---------------------------: " << JacobianFactor[0][0] << std::endl;
    }

// ############################
    if(!d_directions.isSet())
        setDefaultDirections();
    else
        normalizeDirections();
// ###########################
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
// #########################
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


std::pair<double, double> recoverThetaAndPhi(const Eigen::Vector3d& Mu) {
    // Normalizar Mu
    Eigen::Vector3d mu_hat = Mu.normalized();

    // Proyectar sobre el plano XZ y normalizar
    Eigen::Vector3d mu_xz(mu_hat(0), 0, mu_hat(2));
    Eigen::Vector3d mu_xz_hat = mu_xz.normalized();

    // Calcular theta
    double theta = std::acos(mu_hat.dot(mu_xz_hat));
    if (mu_hat(1) > 0) {
        theta = -theta;
    }

    if (theta > 0 && mu_hat(2) < 0) {
        theta = M_PI - theta;
    }

    if (theta < 0 && mu_hat(2) < 0) {
        theta = -(M_PI + theta);
    }

    // Calcular phi
    double phi = std::acos(mu_xz_hat.dot(Eigen::Vector3d(0, 0, 1)));
    if (mu_hat(2) < 0) {
        phi = -std::acos(mu_xz_hat.dot(Eigen::Vector3d(0, 0, -1)));
    }

    if (mu_hat(0) < 0) {
        phi = -phi;
    }

    return {theta, phi};
}

std::pair<double, double> recoverThetaAndPhi_DerivadaZ(const Eigen::Vector3d& Mu) {
    // Normalizar Mu
    Eigen::Vector3d mu_hat = Mu.normalized();
    
    // Calcular theta
    double theta = std::acos(mu_hat.dot(Eigen::Vector3d(0, 0, 1)));
    
    // Proyectar sobre el plano XY y normalizar
    // Verificar esto ------>mu_xy(mu_hat(0), mu_hat(1), 0), ya que en el cuaderno tengo una rotacion de los ejes,
    // por lo que la proyeccion debe ser en xy pero ahora equivalen a xz, o sea esto: mu_xy(mu_hat(1), 0,mu_hat(0))
    Eigen::Vector3d mu_xy(mu_hat(0), mu_hat(1), 0);
    Eigen::Vector3d mu_xy_hat = mu_xy.normalized();
    
    // Calcular phi
    double phi = std::acos(Eigen::Vector3d(0, -1, 0).dot(mu_xy_hat));
    
    if (mu_hat(1) > 0) { // Si el componente Y de Mu es positivo, la rotación es negativa alrededor de X
        theta = -theta;
    }
    
    if (theta < 0) {
        phi = -std::acos(Eigen::Vector3d(0, 1, 0).dot(mu_xy_hat));
    }
    
    if (mu_hat(0) < 0) {
        phi = -phi;
    }
    
    return {theta, phi};
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

    // cout << "Campo_Magnetico_resultado:  " << Campo_Magnetico_resultado.transpose() << endl;


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





    d_constraintIndex.setValue(cIndex);
    const auto& constraintIndex = sofa::helper::getReadAccessor(d_constraintIndex);
    MatrixDeriv& column = *cMatrix.beginEdit();
    const auto& indices = sofa::helper::getReadAccessor(d_indices);
    sofa::Index sizeIndices = indices.size();
    const auto& useDirections = sofa::helper::getReadAccessor(d_useDirections);
    auto directions = sofa::helper::getWriteAccessor(d_directions);
    auto Jacobian = sofa::helper::getWriteAccessor(d_Jacobian);
    auto PosSensor = sofa::helper::getWriteAccessor(d_PosSensor);
    auto mum = sofa::helper::getWriteAccessor(d_mum);
    auto weight = sofa::helper::getReadAccessor(d_weight);

//    std::cout << "datos d epos AAAAAAA: " << x << std::endl;
    auto& data = x;  // Usamos directamente el ReadAccessor    double dBx_drz = 0;
    Eigen::Vector3d B_calculada;  // Variable global o de ámbito extendido
    auto readAccessor = sofa::helper::getReadAccessor(data);  // Obtén un ReadAccessor
    Eigen::Vector3d B_0;  // Variable global o de ámbito extendido
    Eigen::Vector3d B_1;  // Variable global o de ámbito extendido
    Eigen::Vector3d B_2;  // Variable global o de ámbito extendido
    Eigen::Vector3d B_3;  // Variable global o de ámbito extendido
    Eigen::Vector3d B_Phi;  // Variable global o de ámbito extendido
    Eigen::Vector3d B_Phi_devz;  // Variable global o de ámbito extendido
    Eigen::Vector3d dBdR_z;  // Variable global o de ámbito extendido
    Eigen::Vector3d dBdR_x;  // Variable global o de ámbito extendido
    Eigen::Vector3d dBdR_y;  // Variable global o de ámbito extendido
    Eigen::Vector3d dBdTheta;  // Variable global o de ámbito extendido
    Eigen::Vector3d dBdPhi;  // Variable global o de ámbito extendido
    Eigen::Vector3d dBdPhi_devz;  // Variable global o de ámbito extendido

//    Eigen::Vector3d dBdR_y;  // Variable global o de ámbito extendido
    const auto Cambio = 0.01;
    // PosSensor[0] =  sofa::type::Vec<6, double>(0,0,0,0,0,24); //PosSensor se define aqui, debo linkearlo con la escena o algo así, 0,1,2 Para sensor, 3,4,5 Para Iman
    // PosSensor[0] =  sofa::type::Vec<6, double>(2,3,4,5,6,7); //PosSensor se define aqui, debo linkearlo con la escena o algo así
    // PosSensor[1] =  sofa::type::Vec<6, double>(8,9,10,11,12,13); //PosSensor se define aqui, debo linkearlo con la escena o algo así


    for (const auto& coord : readAccessor) {
            // Acceder a todas las componentes de cada Vec<2, double> y mostrar las tres componentes si es posible
            // msg_warning() << "Coord: (" << coord[0] << ", " << coord[1] << ", " << coord[2] << ")";  // Asumiendo que 'Vec<2, double>' tiene tres componentes
            Eigen::Quaterniond MiR(coord[6], coord[3], coord[4], coord[5]);  // (w, x, y, z)
            // std::cout << "MiR de MagnetModel.inl : " << MiR.coeffs().transpose() << "\n";
            Eigen::Matrix3d rotation_matrix = MiR.toRotationMatrix();
            // std::cout << "Matriz de rotación:\n" << rotation_matrix << std::endl;
            // Extraer la última columna
            const double mu_x = rotation_matrix(0, 2);  // Elemento (0, 2)
            const double mu_y = rotation_matrix(1, 2);  // Elemento (1, 2)
            const double mu_z = rotation_matrix(2, 2);  // Elemento (2, 2)
            // std::cout << "mu_x: " << mu_x << ", mu_y: " << mu_y << ", mu_z: " << mu_z << std::endl;
            const double ajuste_x = PosSensor[0][0] ;
            const double ajuste_y = PosSensor[0][1];
            const double ajuste_z = PosSensor[0][2];
            // std::cout << "PosSensor[0] de MagnetModel.inl : " << PosSensor[0] << "\n";
    //        const double ajuste_z = 2.7 - 15; // Este ajuste era para corroborar calculo en c++ con el de python
            // std::cout << "ajuste_x: " << ajuste_x << ", ajuste_y: " << ajuste_y << ", ajuste_z: " << ajuste_z << std::endl;
            // std::cout << "coord - ajuste_x: " << coord[0] - ajuste_x << ", coord - ajuste_y: " << coord[1] - ajuste_y << ", coord - ajuste_z: " << coord[2] - ajuste_z << std::endl;
            
            B_calculada = Calculo_B_Test(coord[0]- ajuste_x,coord[1]- ajuste_y,coord[2] - ajuste_z,mum[0][0],mu_x,mu_y,mu_z);
        //    std::cout << "Campo magnético B_c alculado c++ AAAAAAAAAAAAAA: " << B_calculada.transpose() << std::endl;
            B_0 = Calculo_B_Test(coord[0]- ajuste_x,coord[1]- ajuste_y,coord[2]- ajuste_z,mum[0][0],mu_x,mu_y,mu_z);
            // std::cout << "Campo magnético B_0 calculado c++ AAAAAAAAAAAAAA: " << B_0.transpose() << std::endl;
            // dB/drX
            B_2 = Calculo_B_Test(coord[0] - ajuste_x + Cambio,coord[1]- ajuste_y,coord[2],mum[0][0],mu_x,mu_y,mu_z);
            dBdR_x = (B_2 - B_0)/Cambio;
            // std::cout << "dBdR_x" << dBdR_x.transpose() << std::endl;
            // dB/drY
            B_3 = Calculo_B_Test(coord[0] - ajuste_x,coord[1]- ajuste_y+ Cambio,coord[2]- ajuste_z,mum[0][0],mu_x,mu_y,mu_z);
            dBdR_y = (B_3 - B_0)/Cambio;
            // std::cout << "dBdR_y" << dBdR_y.transpose() << std::endl;
            // dB/drZ
            B_1 = Calculo_B_Test(coord[0] - ajuste_x,coord[1]- ajuste_y,coord[2]- ajuste_z    + Cambio    ,mum[0][0],mu_x,mu_y,mu_z);
            dBdR_z = (B_1 - B_0)/Cambio;
            // std::cout << "dBdR_z" << dBdR_z.transpose() << std::endl;

            
//            ---------------- Angulos :-----------------------------------------------------------

//            Eigen::Vector3d Mu(0.44229157 ,-0.32357449, -0.8364674); // Ejemplo de entrada, Resultado: ThetaRecovered (rad): 2.81208787) PhiRecovered: (rad) -0.4863910100000001
            Eigen::Vector3d Mu(mu_x, mu_y, mu_z );
            auto [theta, phi] = recoverThetaAndPhi(Mu);
            // std::cout << "Theta: " << theta << "\n";
            // std::cout << "Phi: " << phi << "\n";

            // Crear rotaciones alrededor de X y Y
            Eigen::AngleAxisd Rx(theta + Cambio, Eigen::Vector3d::UnitX());
            Eigen::AngleAxisd Ry(phi, Eigen::Vector3d::UnitY());
            Eigen::Quaterniond MiR_2 = Ry * Rx;
            Eigen::Vector3d Mu_theta = MiR_2 * Eigen::Vector3d(0, 0, 1);
            // std::cout << "Mu_theta: " << Mu_theta.transpose() << std::endl;

//            auto [theta_2, phi_2] = recoverThetaAndPhi(Mu_theta);
//            std::cout << "Theta_2: " << theta_2 << "\n";
//            std::cout << "Phi_2: " << phi_2 << "\n";


            // db/dTheta
            B_1 = Calculo_B_Test(coord[0] - ajuste_x ,coord[1]- ajuste_y, coord[2],mum[0][0], Mu_theta[0],Mu_theta[1],Mu_theta[2]);
            dBdTheta= (B_1 - B_0)/ (Cambio);
            // std::cout << "dBdTheta : " << dBdTheta.transpose() << std::endl;

            // db/dPhi
            Eigen::AngleAxisd Rx_Phi(theta, Eigen::Vector3d::UnitX());
            Eigen::AngleAxisd Ry_Phi(phi+ Cambio, Eigen::Vector3d::UnitY());
            Eigen::Quaterniond MiR_Phi = Ry_Phi * Rx_Phi;
            Eigen::Vector3d Mu_Phi = MiR_Phi* Eigen::Vector3d(0, 0, 1);
            B_Phi = Calculo_B_Test(coord[0] - ajuste_x ,coord[1]- ajuste_y, coord[2],mum[0][0], Mu_Phi[0],Mu_Phi[1],Mu_Phi[2]);
            dBdPhi= (B_Phi - B_0)/ (Cambio);
            // std::cout << "dBdPhi : " << dBdPhi.transpose() << std::endl;

//            ---------------- Angulo para eje z :-----------------------------------------------------------
            // Eigen::Vector3d Mu(mu_x, mu_y, mu_z );
            auto [theta_devz, phi_devz] = recoverThetaAndPhi_DerivadaZ(Mu);
            // std::cout << "Theta_devz: " << theta_devz << "\n";
            // std::cout << "Phi_devz: " << phi_devz << "\n";     

            // Crear rotaciones alrededor de X y Y
            Eigen::AngleAxisd Rx_devz(theta_devz , Eigen::Vector3d::UnitX());
            Eigen::AngleAxisd Rz_devz(phi_devz + Cambio, Eigen::Vector3d::UnitZ());
            Eigen::Quaterniond MiR_devz = Rz_devz * Rx_devz;
            Eigen::Vector3d Mu_devz = MiR_devz * Eigen::Vector3d(0, 0, 1);

            B_Phi_devz = Calculo_B_Test(coord[0] - ajuste_x ,coord[1]- ajuste_y, coord[2],mum[0][0], Mu_devz[0],Mu_devz[1],Mu_devz[2]);
            // std::cout << "mum de MagnetModel.inl : " << mum[0][0] << "\n";     
            dBdPhi_devz= (B_Phi_devz - B_0)/ (Cambio);
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

    double dBx_dTheta = dBdTheta[0];
    double dBy_dTheta = dBdTheta[1];
    double dBz_dTheta = dBdTheta[2];

    double dBx_dPhi = dBdPhi[0];
    double dBy_dPhi = dBdPhi[1];
    double dBz_dPhi = dBdPhi[2];

    double dBx_dPhi_devz = dBdPhi_devz[0];
    double dBy_dPhi_devz = dBdPhi_devz[1];
    double dBz_dPhi_devz = dBdPhi_devz[2];

//    Jacobian.push_back(VecDeriv(0.0, 0.0, dBx_drz, 0.0, 0.0, 0.0));


//    Jacobian[0] = sofa::type::Vec<6, double>(17114, 0.0, 0.0, 0.0, 0.0, 0.0);  // Asignar un nuevo valor
//    Jacobian[1] = sofa::type::Vec<6, double>(0.0, 17114, 0.0, 0.0, 0.0, 0.0);  // Asignar un nuevo valor
//    Jacobian[2] = sofa::type::Vec<6, double>(-126.7, -126.7, -33978, 0.0, 0.0, 0.0);  // Asignar un nuevo valor

        // double f = 100.0; // Prueba subiendo este valor a 1000 o 10000 si sigue lento
        // double f = 1.0; // Prueba subiendo este valor a 1000 o 10000 si sigue lento
        double f = d_JacobianFactor.getValue()[0][0];
        // std::cout << "Actual JacobianFactor " << f << std::endl;
        Jacobian[0] = sofa::type::Vec<6, double>(dBx_drx/f, dBx_dry/f, dBx_drz/f, dBx_dTheta/f, dBx_dPhi/f, dBx_dPhi_devz/f);
        Jacobian[1] = sofa::type::Vec<6, double>(dBy_drx/f, dBy_dry/f, dBy_drz/f, dBy_dTheta/f, dBy_dPhi/f, dBy_dPhi_devz/f);
        Jacobian[2] = sofa::type::Vec<6, double>(dBz_drx/f, dBz_dry/f, dBz_drz/f, dBz_dTheta/f, dBz_dPhi/f, dBz_dPhi_devz/f);


    unsigned int index = 0;

    for (unsigned j = 0; j < 3; j++) {
        MatrixDerivRowIterator rowIterator = column.writeLine(constraintIndex + index);
    //    std::cout << "datos de directions: " << directions[j] << std::endl;

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
    // d_PosSensor.setValue(directions);
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

