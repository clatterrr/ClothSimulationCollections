#pragma once
#include <Eigen\Core>
#include <Eigen\Dense>
using namespace Eigen;

#define SPARSE_DF

const int Nx = 2;
const int Ny = 2;
const int node_num = Nx * Ny;
VectorXf forces(node_num * 3);

class Condition
{
public:
	float alpha_shear;
	float k_shear;
	float damping_shear;

	float alpha_stretch;
	float k_stretch;
	float damping_stretch;
	Matrix3f identity;

	Condition()
	{

		alpha_shear = 1;
		k_shear = 1;
		damping_shear = 2;

		alpha_stretch = 1;
		k_stretch = 1;
		damping_stretch = 2;

		identity << 1, 0, 0, 0, 1, 0, 0, 0, 1;
	}



public:


	void ComputeShearForces(Vector3f wu, Vector3f wv, Vector3f dwu, Vector3f dwv, Vector3f p0, Vector3f p1, Vector3f p2, Vector3f v0, Vector3f v1, Vector3f v2, int idx0, int idx1, int idx2)//, MatrixXf dfdx)
	{

		float ShearCondition = alpha_shear * wu.dot(wv);


		Vector3f dcdx0 = alpha_shear * (dwu(0) * wv + dwv(0) * wu);
		Vector3f dcdx1 = alpha_shear * (dwu(1) * wv + dwv(1) * wu);
		Vector3f dcdx2 = alpha_shear * (dwu(2) * wv + dwv(2) * wu);

		Matrix3f d2dx0x0 = alpha_shear * (dwu(0) * dwv(0) + dwu(0) * dwu(0)) * identity;
		Matrix3f d2dx0x1 = alpha_shear * (dwu(0) * dwv(1) + dwu(1) * dwu(0)) * identity;
		Matrix3f d2dx0x2 = alpha_shear * (dwu(0) * dwv(2) + dwu(2) * dwu(0)) * identity;

		Matrix3f d2dx1x0 = alpha_shear * (dwu(1) * dwv(0) + dwu(0) * dwu(1)) * identity;
		Matrix3f d2dx1x1 = alpha_shear * (dwu(1) * dwv(1) + dwu(1) * dwu(1)) * identity;
		Matrix3f d2dx1x2 = alpha_shear * (dwu(1) * dwv(2) + dwu(2) * dwu(1)) * identity;

		Matrix3f d2dx2x0 = alpha_shear * (dwu(2) * dwv(1) + dwu(1) * dwu(2)) * identity;
		Matrix3f d2dx2x1 = alpha_shear * (dwu(2) * dwv(1) + dwu(1) * dwu(2)) * identity;
		Matrix3f d2dx2x2 = alpha_shear * (dwu(2) * dwv(2) + dwu(2) * dwu(2)) * identity;

		// ¼ÆËãÁ¦

		forces.segment<3>(idx0 * 3) += -k_shear * ShearCondition * dcdx0;
		forces.segment<3>(idx1 * 3) += -k_shear * ShearCondition * dcdx1;
		forces.segment<3>(idx2 * 3) += -k_shear * ShearCondition * dcdx2;
	}
	void ComputeStretchForces(Vector3f wu, Vector3f wv, Vector3f dwu, Vector3f dwv, Vector3f p0, Vector3f p1, Vector3f p2, Vector3f v0, Vector3f v1, Vector3f v2, int idx0, int idx1, int idx2)//, MatrixXf dfdx)
	{

		float wuNorm = wu.norm();
		float wvNorm = wv.norm();

		float cu = alpha_stretch * (wuNorm - 1);
		float cv = alpha_stretch * (wvNorm - 1);

		Vector3f dcudx0 = alpha_stretch * dwu(0) * wu / wuNorm;
		Vector3f dcudx1 = alpha_stretch * dwu(1) * wu / wuNorm;
		Vector3f dcudx2 = alpha_stretch * dwu(2) * wu / wuNorm;

		Vector3f dcvdx0 = alpha_stretch * dwv(0) * wv / wvNorm;
		Vector3f dcvdx1 = alpha_stretch * dwv(1) * wv / wvNorm;
		Vector3f dcvdx2 = alpha_stretch * dwv(2) * wv / wvNorm;

		forces.segment<3>(idx0 * 3) += -k_stretch * (cu * dcudx0 + cv * dcvdx0);
		forces.segment<3>(idx1 * 3) += -k_stretch * (cu * dcudx1 + cv * dcvdx1);
		forces.segment<3>(idx2 * 3) += -k_stretch * (cu * dcudx2 + cv * dcvdx2);

	}
};