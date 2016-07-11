#ifndef ELTENSORS_H
#define ELTENSORS_H
#include "ELInitialization.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
class ELTensors {

private:
	//4-tensor of Cell X
	//for example CellX[0][0][1]=X(x,y,z+1) and CellX[0][0][2]=X(x,y,z-1)	
	std::vector<std::vector<std::vector<std::vector<double> > > > CellX;

	//mesh lengths
	double dx;
	double dy;
	double dz;
	//volume fraction
	double Volfn;
	ELInitialization Eli;
	std::vector<double> Ahat;

public:
	std::ofstream& myfile;

	/// Constructor.
	ELTensors(std::vector<std::vector<std::vector<std::vector<double>>>> xCellX,
		double xdx,
		double xdy,
		double xdz,
		double xVolfn,
		ELInitialization xEli,
		std::vector<double> xAhat,
		std::ofstream& xmyfile
		) : CellX(xCellX), dx(xdx), dy(xdy), dz(xdz), Volfn(xVolfn), Eli(xEli), Ahat(xAhat),myfile(xmyfile) { };


	/// Copy constructor.



	//Member Functions

	//ELTensors

	//Embedded metric h_ij
	const std::vector<double> hmatrix();

	//Diffusion tensor product with diffusion direction q: qDq
	const double qDq(std::vector<double> qk);


	//Gamma^i_{j,k}=1/2 h^(il)[partial_{j}h_{lk}+ partial_{k}h_{jl}-partial_{l}h_{jk}]
	const double CSIComp(double i, double j, double k) const;

	//entries of gamma
	const std::vector<std::vector<double>> IMmatrix();

	//determinant gamma
	const double detgamma();

	//Induced metric gamma_nu,mu
	const std::vector<std::vector<double>> IIMComp();


	//Tensor derivatives

	//Derivative of embedding map X
	const std::vector<std::vector<double>> DerivX();

	//Second derivative of embedding map X
	const std::vector<std::vector<std::vector<double>>> DDerivX();

	//Derivative of induced metric gamma
	const std::vector<std::vector<std::vector<double>>> Derivgamma();

	//Deriv determinant gamma
	const double Derivdetgamma(double direction);


	//Derivative of inverse induced metric gamma_nu,mu
	const std::vector<std::vector<double>> DerivIIMComp();

	//Product of Derivative of Diffusion tensor D with kth diffusion direction q_k: (q_k^T) partialD (q_k)
	const double qpartialDq(double direction, std::vector<double> qk);


	//Euler Lagrange equations and Volume fraction

	///Our EL eqns have three terms: 1)The sum involving the diffusion tensors , 
	const double term1Diff(double direction);

	/// 2)the partial derivatives of induced metric gamma with embedding map X and
	const double term2IX(double direction);

	///3) the term involving the Christoffel symbols.
	const double term3CS(double direction);

	///We will compute them separately and then add them for the EL scheme
	const double ELequation(double direction);

	//The iteration rule for the volume fraction
	const double VolfraIter();


};







#endif

