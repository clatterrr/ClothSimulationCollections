#pragma once
#include <Eigen\Core>
#include <Eigen\Dense>
using namespace Eigen;



const int Nx = 3;
const int Ny = 3;
const int node_num = Nx * Ny;


Matrix<float, node_num * 3, node_num * 3> dfdx;
Matrix<float, node_num * 3, node_num * 3> dfdv;
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


		df0dx0 = -k_shear * (dcdx0 * dcdx0.transpose() + ShearCondition * d2dx0x0);
		df0dx1 = -k_shear * (dcdx0 * dcdx1.transpose() + ShearCondition * d2dx0x1);
		df0dx2 = -k_shear * (dcdx0 * dcdx2.transpose() + ShearCondition * d2dx0x2);


		df1dx0 = -k_shear * (dcdx1 * dcdx0.transpose() + ShearCondition * d2dx1x0);
		df1dx1 = -k_shear * (dcdx1 * dcdx1.transpose() + ShearCondition * d2dx1x1);
		df1dx2 = -k_shear * (dcdx1 * dcdx2.transpose() + ShearCondition * d2dx1x2);

		Matrix3f df2dx0 = -k_shear * (dcdx2 * dcdx0.transpose() + ShearCondition * d2dx2x0);
		Matrix3f df2dx1 = -k_shear * (dcdx2 * dcdx1.transpose() + ShearCondition * d2dx2x1);
		Matrix3f df2dx2 = -k_shear * (dcdx2 * dcdx2.transpose() + ShearCondition * d2dx2x2);

		dfdx.block<3, 3>(idx0 * 3, idx0 * 3) += df0dx0;
		dfdx.block<3, 3>(idx0 * 3, idx1 * 3) += df0dx1;
		dfdx.block<3, 3>(idx0 * 3, idx2 * 3) += df0dx2;

		dfdx.block<3, 3>(idx1 * 3, idx0 * 3) += df1dx0;
		dfdx.block<3, 3>(idx1 * 3, idx1 * 3) += df1dx1;
		dfdx.block<3, 3>(idx1 * 3, idx2 * 3) += df1dx2;

		dfdx.block<3, 3>(idx2 * 3, idx0 * 3) += df2dx0;
		dfdx.block<3, 3>(idx2 * 3, idx1 * 3) += df2dx1;
		dfdx.block<3, 3>(idx2 * 3, idx2 * 3) += df2dx2;



		df0dv0 = -damping_shear * (dcdx0 * dcdx0.transpose());
		df0dv1 = -damping_shear * (dcdx0 * dcdx1.transpose());
		df0dv2 = -damping_shear * (dcdx0 * dcdx2.transpose());

		df1dv0 = -damping_shear * (dcdx1 * dcdx0.transpose());
		df1dv1 = -damping_shear * (dcdx1 * dcdx1.transpose());
		df1dv2 = -damping_shear * (dcdx1 * dcdx2.transpose());

		df2dv0 = -damping_shear * (dcdx2 * dcdx0.transpose());
		df2dv1 = -damping_shear * (dcdx2 * dcdx1.transpose());
		df2dv2 = -damping_shear * (dcdx2 * dcdx2.transpose());

		dfdv.block<3, 3>(idx0 * 3, idx0 * 3) += df0dv0;
		dfdv.block<3, 3>(idx0 * 3, idx1 * 3) += df0dv1;
		dfdv.block<3, 3>(idx0 * 3, idx2 * 3) += df0dv2;

		dfdv.block<3, 3>(idx1 * 3, idx0 * 3) += df1dv0;
		dfdv.block<3, 3>(idx1 * 3, idx1 * 3) += df1dv1;
		dfdv.block<3, 3>(idx1 * 3, idx2 * 3) += df1dv2;

		dfdv.block<3, 3>(idx2 * 3, idx0 * 3) += df2dv0;
		dfdv.block<3, 3>(idx2 * 3, idx1 * 3) += df2dv1;
		dfdv.block<3, 3>(idx2 * 3, idx2 * 3) += df2dv2;

		df0dx0 = -damping_shear * (d2dx0x0 * dcdt);
		df0dx1 = -damping_shear * (d2dx0x1 * dcdt);
		df0dx2 = -damping_shear * (d2dx0x2 * dcdt);


		df1dx0 = -damping_shear * (d2dx1x0 * dcdt);
		df1dx1 = -damping_shear * (d2dx1x1 * dcdt);
		df1dx2 = -damping_shear * (d2dx1x2 * dcdt);

		df2dx0 = -damping_shear * (d2dx2x0 * dcdt);
		df2dx1 = -damping_shear * (d2dx2x1 * dcdt);
		df2dx2 = -damping_shear * (d2dx2x2 * dcdt);

		dfdx.block<3, 3>(idx0 * 3, idx0 * 3) += df0dx0;
		dfdx.block<3, 3>(idx0 * 3, idx1 * 3) += df0dx1;
		dfdx.block<3, 3>(idx0 * 3, idx2 * 3) += df0dx2;

		dfdx.block<3, 3>(idx1 * 3, idx0 * 3) += df1dx0;
		dfdx.block<3, 3>(idx1 * 3, idx1 * 3) += df1dx1;
		dfdx.block<3, 3>(idx1 * 3, idx2 * 3) += df1dx2;

		dfdx.block<3, 3>(idx2 * 3, idx0 * 3) += df2dx0;
		dfdx.block<3, 3>(idx2 * 3, idx1 * 3) += df2dx1;
		dfdx.block<3, 3>(idx2 * 3, idx2 * 3) += df2dx2;


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



		Matrix3f wuMatrix = (Matrix3f::Identity() - wu * wu.transpose() / (wuNorm * wuNorm));
		Matrix3f wvMatrix = (Matrix3f::Identity() - wv * wv.transpose() / (wvNorm * wvNorm));

		Matrix3f d2cudx0x0 = (alpha_stretch / wuNorm) * dwu(0) * dwu(0) * wuMatrix;
		Matrix3f d2cudx0x1 = (alpha_stretch / wuNorm) * dwu(0) * dwu(1) * wuMatrix;
		Matrix3f d2cudx0x2 = (alpha_stretch / wuNorm) * dwu(0) * dwu(2) * wuMatrix;

		Matrix3f d2cudx1x0 = (alpha_stretch / wuNorm) * dwu(1) * dwu(0) * wuMatrix;
		Matrix3f d2cudx1x1 = (alpha_stretch / wuNorm) * dwu(1) * dwu(1) * wuMatrix;
		Matrix3f d2cudx1x2 = (alpha_stretch / wuNorm) * dwu(1) * dwu(2) * wuMatrix;

		Matrix3f d2cudx2x0 = (alpha_stretch / wuNorm) * dwu(2) * dwu(0) * wuMatrix;
		Matrix3f d2cudx2x1 = (alpha_stretch / wuNorm) * dwu(2) * dwu(1) * wuMatrix;
		Matrix3f d2cudx2x2 = (alpha_stretch / wuNorm) * dwu(2) * dwu(2) * wuMatrix;

		Matrix3f d2cvdx0x0 = (alpha_stretch / wvNorm) * dwv(0) * dwv(0) * wvMatrix;
		Matrix3f d2cvdx0x1 = (alpha_stretch / wvNorm) * dwv(0) * dwv(1) * wvMatrix;
		Matrix3f d2cvdx0x2 = (alpha_stretch / wvNorm) * dwv(0) * dwv(2) * wvMatrix;

		Matrix3f d2cvdx1x0 = (alpha_stretch / wvNorm) * dwv(1) * dwv(0) * wvMatrix;
		Matrix3f d2cvdx1x1 = (alpha_stretch / wvNorm) * dwv(1) * dwv(1) * wvMatrix;
		Matrix3f d2cvdx1x2 = (alpha_stretch / wvNorm) * dwv(1) * dwv(2) * wvMatrix;

		Matrix3f d2cvdx2x0 = (alpha_stretch / wvNorm) * dwv(2) * dwv(0) * wvMatrix;
		Matrix3f d2cvdx2x1 = (alpha_stretch / wvNorm) * dwv(2) * dwv(1) * wvMatrix;
		Matrix3f d2cvdx2x2 = (alpha_stretch / wvNorm) * dwv(2) * dwv(2) * wvMatrix;

		df0dx0 = -k_stretch * (dcudx0 * dcudx0.transpose() + cu * d2cudx0x0 + dcvdx0 * dcvdx0.transpose() + cv * d2cvdx0x0);
		df0dx1 = -k_stretch * (dcudx0 * dcudx1.transpose() + cu * d2cudx0x1 + dcvdx0 * dcvdx1.transpose() + cv * d2cvdx0x1);
		df0dx2 = -k_stretch * (dcudx0 * dcudx2.transpose() + cu * d2cudx0x2 + dcvdx0 * dcvdx2.transpose() + cv * d2cvdx0x2);

		df1dx0 = -k_stretch * (dcudx1 * dcudx0.transpose() + cu * d2cudx1x0 + dcvdx1 * dcvdx0.transpose() + cv * d2cvdx1x0);
		df1dx1 = -k_stretch * (dcudx1 * dcudx1.transpose() + cu * d2cudx1x1 + dcvdx1 * dcvdx1.transpose() + cv * d2cvdx1x1);
		df1dx2 = -k_stretch * (dcudx1 * dcudx2.transpose() + cu * d2cudx1x2 + dcvdx1 * dcvdx2.transpose() + cv * d2cvdx1x2);

		df2dx0 = -k_stretch * (dcudx2 * dcudx0.transpose() + cu * d2cudx2x0 + dcvdx2 * dcvdx0.transpose() + cv * d2cvdx2x0);
		df2dx1 = -k_stretch * (dcudx2 * dcudx1.transpose() + cu * d2cudx2x1 + dcvdx2 * dcvdx1.transpose() + cv * d2cvdx2x1);
		df2dx2 = -k_stretch * (dcudx2 * dcudx2.transpose() + cu * d2cudx2x2 + dcvdx2 * dcvdx2.transpose() + cv * d2cvdx2x2);

		dfdx.block<3, 3>(idx0 * 3, idx0 * 3) += df0dx0;
		dfdx.block<3, 3>(idx0 * 3, idx1 * 3) += df0dx1;
		dfdx.block<3, 3>(idx0 * 3, idx2 * 3) += df0dx2;

		dfdx.block<3, 3>(idx1 * 3, idx0 * 3) += df1dx0;
		dfdx.block<3, 3>(idx1 * 3, idx1 * 3) += df1dx1;
		dfdx.block<3, 3>(idx1 * 3, idx2 * 3) += df1dx2;

		dfdx.block<3, 3>(idx2 * 3, idx0 * 3) += df2dx0;
		dfdx.block<3, 3>(idx2 * 3, idx1 * 3) += df2dx1;
		dfdx.block<3, 3>(idx2 * 3, idx2 * 3) += df2dx2;


		// 速度阻尼项
		float dcudt = dcudx0.dot(v0) + dcudx1.dot(v1) + dcudx2.dot(v2);
		float dcvdt = dcvdx0.dot(v0) + dcvdx1.dot(v1) + dcvdx2.dot(v2);

		forces.segment<3>(idx0 * 3) += -damping_stretch * (dcudt * dcudx0 + dcvdt * dcvdx0);
		forces.segment<3>(idx1 * 3) += -damping_stretch * (dcudt * dcudx1 + dcvdt * dcvdx1);
		forces.segment<3>(idx2 * 3) += -damping_stretch * (dcudt * dcudx2 + dcvdt * dcvdx2);

		df0dv0 = -damping_stretch * (dcudx0 * dcudx0.transpose() + dcvdx0 * dcvdx0.transpose());
		df0dv1 = -damping_stretch * (dcudx0 * dcudx1.transpose() + dcvdx0 * dcvdx1.transpose());
		df0dv2 = -damping_stretch * (dcudx0 * dcudx2.transpose() + dcvdx0 * dcvdx2.transpose());

		df1dv0 = -damping_stretch * (dcudx1 * dcudx0.transpose() + dcvdx1 * dcvdx0.transpose());
		df1dv1 = -damping_stretch * (dcudx1 * dcudx1.transpose() + dcvdx1 * dcvdx1.transpose());
		df1dv2 = -damping_stretch * (dcudx1 * dcudx2.transpose() + dcvdx1 * dcvdx2.transpose());

		df2dv0 = -damping_stretch * (dcudx2 * dcudx0.transpose() + dcvdx2 * dcvdx0.transpose());
		df2dv1 = -damping_stretch * (dcudx2 * dcudx1.transpose() + dcvdx2 * dcvdx1.transpose());
		df2dv2 = -damping_stretch * (dcudx2 * dcudx2.transpose() + dcvdx2 * dcvdx2.transpose());

		dfdv.block<3, 3>(idx0 * 3, idx0 * 3) += df0dv0;
		dfdv.block<3, 3>(idx0 * 3, idx1 * 3) += df0dv1;
		dfdv.block<3, 3>(idx0 * 3, idx2 * 3) += df0dv2;

		dfdv.block<3, 3>(idx1 * 3, idx0 * 3) += df1dv0;
		dfdv.block<3, 3>(idx1 * 3, idx1 * 3) += df1dv1;
		dfdv.block<3, 3>(idx1 * 3, idx2 * 3) += df1dv2;

		dfdv.block<3, 3>(idx2 * 3, idx0 * 3) += df2dv0;
		dfdv.block<3, 3>(idx2 * 3, idx1 * 3) += df2dv1;
		dfdv.block<3, 3>(idx2 * 3, idx2 * 3) += df2dv2;


		df0dx0 = -damping_stretch * (d2cudx0x0 * dcudt + d2cvdx0x0 * dcvdt);
		df0dx1 = -damping_stretch * (d2cudx0x1 * dcudt + d2cvdx0x1 * dcvdt);
		df0dx2 = -damping_stretch * (d2cudx0x2 * dcudt + d2cvdx0x2 * dcvdt);


		df1dx0 = -damping_stretch * (d2cudx1x0 * dcudt + d2cvdx1x0 * dcvdt);
		df1dx1 = -damping_stretch * (d2cudx1x1 * dcudt + d2cvdx1x1 * dcvdt);
		df1dx2 = -damping_stretch * (d2cudx1x2 * dcudt + d2cvdx1x2 * dcvdt);

		df2dx0 = -damping_stretch * (d2cudx2x0 * dcudt + d2cvdx2x0 * dcvdt);
		df2dx1 = -damping_stretch * (d2cudx2x1 * dcudt + d2cvdx2x1 * dcvdt);
		df2dx2 = -damping_stretch * (d2cudx2x2 * dcudt + d2cvdx2x2 * dcvdt);

		dfdx.block<3, 3>(idx0 * 3, idx0 * 3) += df0dx0;
		dfdx.block<3, 3>(idx0 * 3, idx1 * 3) += df0dx1;
		dfdx.block<3, 3>(idx0 * 3, idx2 * 3) += df0dx2;

		dfdx.block<3, 3>(idx1 * 3, idx0 * 3) += df1dx0;
		dfdx.block<3, 3>(idx1 * 3, idx1 * 3) += df1dx1;
		dfdx.block<3, 3>(idx1 * 3, idx2 * 3) += df1dx2;

		dfdx.block<3, 3>(idx2 * 3, idx0 * 3) += df2dx0;
		dfdx.block<3, 3>(idx2 * 3, idx1 * 3) += df2dx1;
		dfdx.block<3, 3>(idx2 * 3, idx2 * 3) += df2dx2;


	}
};