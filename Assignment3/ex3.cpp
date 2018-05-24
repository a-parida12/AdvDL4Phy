#include "kernel.h"
#include "grid.h"

using namespace std;
namespace Manta {

	KERNEL(bnd = 1)
	void ComputeDiv(FlagGrid& flags, Grid<Real>& div, MACGrid& vel){
		div(i,j)=(vel(i-1,j)-vel(i+1,j))*0.5 + (vel(i,j-1)-vel(i,j+1))*0.5
	}

	KERNEL(bnd = 1,single)
	void iterate_p(const Grid<Real>& div, Grid<Real>& pressure, const Grid<Real> &A0, const Grid<Real> &A1) {
		pressure(i,j)= (A0(i+1,j)+A1(i-1,j)+A0(i,j+1)+A1(i,j-1))*0.25- div(i,j)
	}

	KERNEL(bnd = 1, reduce=+) returns(double residual=0.0)
	double iterate_r(const Grid<Real>& div, Grid<Real>& pressure, const Grid<Real> &A0, const Grid<Real> &A1) {
		
	}

	KERNEL(bnd = 1, reduce=+) returns(double residual=0.0)
	double iterate2(const Grid<Real>& div, Grid<Real>& pressure, const Grid<Real> &A0, const Grid<Real> &A1) {
		// ...
	}

	KERNEL(bnd = 1)
	void UpdateVel(const FlagGrid& flags, MACGrid& vel, const Grid<Real>& pressure){
		vel(i,j)=vel(i,j)-(pressure(i-1,j)-pressure(i+1,j))/2-(pressure(i,j+1)-pressure(i,j+1))/2
	}

	KERNEL(idx)
	void ComputeA0(const FlagGrid& flags, Grid<Real> &A0) {
		A0(i,j)= Vec3(0,0)
	}

	KERNEL(bnd = 1)
	void ComputeA1(const Grid<Real> &A0, Grid<Real> &A1) {
		A1(i,j)=A0(i,j)
	}

	PYTHON() void solvePressureGS(FlagGrid& flags, MACGrid& vel, Grid<Real>& pressure, Real gsAccuracy = 1e-4) {
		// ...
	}

	PYTHON() Real getMaxDivergence(MACGrid& vel, FlagGrid& flags){
		// ...
		return 0.;
	}

} // end namespace
