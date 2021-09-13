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

	Vector3f v10, v20, v23, v13, v12;
	Vector3f e10, e20, e23, e13, e12;
	float c00, c01, c02, c13, c11, c12;
	Vector3f n0, n1;
	Vector3f b00, b01, b02, b13, b11, b12;
	float d00, d01, d02, d13, d11, d12;
	Matrix3f identity;

	// Bend
	float alpha_bend;
	float k_bend;
	float damping_bend;

	Vector3f dThetadx0, dThetadx1, dThetadx2, dThetadx3;

	float theta, sinTheta, cosTheta, dThetadt;

	Matrix3f dn0dx0, dn0dx1, dn0dx2, dn0dx3, dn1dx0, dn1dx1, dn1dx2, dn1dx3;
	// derivatives of the cosines:
	Vector3f dc01dx0, dc01dx1, dc01dx2, dc01dx3, dc02dx0, dc02dx1, dc02dx2, dc02dx3;
	Vector3f dc11dx0, dc11dx1, dc11dx2, dc11dx3, dc12dx0, dc12dx1, dc12dx2, dc12dx3;
	// derivatives of the perpendicular distances:
	Vector3f dd00dx0, dd00dx1, dd00dx2, dd00dx3;
	Vector3f dd01dx0, dd01dx1, dd01dx2, dd01dx3;
	Vector3f dd02dx0, dd02dx1, dd02dx2, dd02dx3;
	Vector3f dd11dx0, dd11dx1, dd11dx2, dd11dx3;
	Vector3f dd12dx0, dd12dx1, dd12dx2, dd12dx3;
	Vector3f dd13dx0, dd13dx1, dd13dx2, dd13dx3;
	// second derivatives of theta with respect to the different vertex positions:
	Matrix3f d2Thetadx0dx0, d2Thetadx0dx1, d2Thetadx0dx2, d2Thetadx0dx3;
	Matrix3f d2Thetadx1dx0, d2Thetadx1dx1, d2Thetadx1dx2, d2Thetadx1dx3;
	Matrix3f d2Thetadx2dx0, d2Thetadx2dx1, d2Thetadx2dx2, d2Thetadx2dx3;
	Matrix3f d2Thetadx3dx0, d2Thetadx3dx1, d2Thetadx3dx2, d2Thetadx3dx3;

	// 计算力用到的
	Matrix3f df0dx0, df0dx1, df0dx2, df0dx3;
	Matrix3f df1dx0, df1dx1, df1dx2, df1dx3;
	Matrix3f df2dx0, df2dx1, df2dx2, df2dx3;
	Matrix3f df3dx0, df3dx1, df3dx2, df3dx3;

	Matrix3f df0dv0, df0dv1, df0dv2, df0dv3;
	Matrix3f df1dv0, df1dv1, df1dv2, df1dv3;
	Matrix3f df2dv0, df2dv1, df2dv2, df2dv3;
	Matrix3f df3dv0, df3dv1, df3dv2, df3dv3;

	Condition()
	{

		alpha_shear = 1;
		k_shear = 1;
		damping_shear = 2;

		alpha_stretch = 1;
		k_stretch = 1;
		damping_stretch = 2;

		alpha_bend = 1;
		k_bend = 1;
		damping_bend = 2;

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

		// 计算力

		forces.segment<3>(idx0 * 3) += -k_shear * ShearCondition * dcdx0;
		forces.segment<3>(idx1 * 3) += -k_shear * ShearCondition * dcdx1;
		forces.segment<3>(idx2 * 3) += -k_shear * ShearCondition * dcdx2;


		float dcdt = dcdx0.dot(v0) + dcdx1.dot(v1) + dcdx2.dot(v2);

		forces.segment<3>(idx0 * 3) += -damping_shear * dcdt * dcdx0;
		forces.segment<3>(idx1 * 3) += -damping_shear * dcdt * dcdx1;
		forces.segment<3>(idx2 * 3) += -damping_shear * dcdt * dcdx2;




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

		// 速度阻尼项
		float dcudt = dcudx0.dot(v0) + dcudx1.dot(v1) + dcudx2.dot(v2);
		float dcvdt = dcvdx0.dot(v0) + dcvdx1.dot(v1) + dcvdx2.dot(v2);

		forces.segment<3>(idx0 * 3) += -damping_stretch * (dcudt * dcudx0 + dcvdt * dcvdx0);
		forces.segment<3>(idx1 * 3) += -damping_stretch * (dcudt * dcudx1 + dcvdt * dcvdx1);
		forces.segment<3>(idx2 * 3) += -damping_stretch * (dcudt * dcudx2 + dcvdt * dcvdx2);

	}




};