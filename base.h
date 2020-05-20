#pragma once
#ifndef BASE_H
#define BASE_H

#include <BaseSDC.h>
#include <BaseMatrix.h>
#include <SiftGPU\SiftGPU.h>
#include <BaseMath.h>
#include "Vocabulary.h"
using namespace basetk;

typedef DBoW3::Vocabulary SIFTVocabulary;

typedef struct tagSIFTInfo
{

	vector<SiftGPU::SiftKeypoint> Keys;
	vector<float> UR;
	vector<unsigned char> des;
	vector<long int> Point_index;
	int KeyNum;
	tagSIFTInfo()
	{
		ZeroStruct(*this, tagSIFTInfo);

	}
}SIFTInfo;

typedef struct tagUV
{
	int Frame_id;
	int Des_id;
	double UV[3];
	int ValidFlag;
	double BackProjectUVError[3];
	double MeanBP_UVError;
	tagUV()
	{
		ZeroStruct(*this, tagUV);

	}
}UV;

typedef struct tagPointXYZ
{
	long int Point_id;
	double XYZ[3];
	double VarXYZ[9];
	int Recent_FI;//要更新
	int n_Frame;
	vector<UV> uvlist;

	tagPointXYZ()
	{
		ZeroStruct(*this, tagPointXYZ);

	}
}PointXYZ;


typedef struct tagPInfo
{
	double MMS_Lat[3] = { 0.0 };
	double MMS_Dis[3] = { 0.0 };
	double MMS_Ang[3] = { 0.0 };
	double MMS_Dep[3] = { 0.0 };
	double Res = 0.0;
	double Post_sigma2=0.0;
	double Dop;
	tagPInfo()
	{
		ZeroStruct(*this, tagPInfo);
	}
}PInfo;


typedef struct tagGroundTruth
{
	GPSTIME gt;
	double Pos[3];
	double VarPos[9];
	double PPos[9];
	double Att[3];
	double VarAtt[9];//对应的失准角方差
	double PAtt[9];
	double Rwc[9];
	int AmbState;
	tagGroundTruth()
	{
		ZeroStruct(*this, tagGroundTruth);

	}
}GroundTruth;

typedef struct tagDes_t
{
	vector<uchar> des;
	int GroupID;		//聚类中心ID

	tagDes_t()
	{
		ZeroStruct(*this, tagDes_t);
	}
	tagDes_t(vector<uchar> &p, int id)
	{
		des = p;
		GroupID = id;
	}
}Des_t;


//普通的通用的最小二乘
class CLSEstimator
{

public:

	CLSEstimator();
	virtual ~CLSEstimator();

	bool initLSQ(int nState, int nMeas, bool bPDiagonal); ///< init LSQ
	bool LSQ(bool bPostEst = false); ///< Least-Square Estimate X=(BTPB)^-1*BTPL

	// 输入
	Mat m_B;      ///< 观测矩阵
	Mat m_P;      ///< 观测权阵,不是协防差阵
	Mat m_L;      ///< 观测值

	// 输出
	Mat m_StateX; ///< 估计的状态
	Mat m_StateP; ///< 估计的状态协方差
	Mat m_V;      ///< 改正数 V = BX - W

	double m_Sigma;

private:

	bool m_bInit;      ///< 是否初始化LSQ
	bool m_bPDiagonal; ///< 观测权阵P是否为对角线
	int  m_nState;     ///< 状态数
	int  m_nMeas;      ///< 观测数
	Mat  m_W;          ///< BTPL
	double m_SigmaPre; ///< 当前sigma不能用时,用前一个代入

	void FreeLSQ();
	void ClearLSQ();
};




class CLSEstGeneral
{

public:

	CLSEstGeneral();
	virtual ~CLSEstGeneral();

	bool initLSQ(int nStateT, int DimsT, int nStateX, int DimsX, int nMeas, int DimsMeas);
	bool AddInfo(int id, Mat A, Mat B, Mat L, Mat P, Mat PriP_T, Mat PriLT);
	bool AddPriInfo(Mat PriP_X, Mat PriLx);
	bool LSQ();

	//输入
	Mat m_N11;
	Mat m_N11i;
	Mat m_N22;
	Mat m_N12;
	Mat m_W1;
	Mat m_W2;
	Mat m_N12TN11iN12;
	Mat m_N12TN11iW1;

	//输出
	Mat m_StateT;
	Mat m_StateX;
	Mat m_StateTP;
	Mat m_StateXP;

private:
	int m_nStateT;
	int m_DimsT;
	int m_nStateX;
	int m_DimsX;
	int m_nMeas;
	int m_DimsMeas;

	void FreeLSQ();
	void ClearLSQ();
};

class CKMeans
{
public:
	CKMeans();
	virtual ~CKMeans();

	int m_k;

	typedef vector<tagDes_t> VecDes_t;

	VecDes_t m_desDB;				//待聚类的特征子descriptor databese
	vector<VecDes_t> m_grp_desDB;	//K类，每一类存储若干特征子
	vector<vector<uchar>> m_center; //聚类中心

	inline void setK(int k_)
	{
		m_k = k_;
		m_grp_desDB.resize(m_k);
	}
	
	bool SetInputDesDatabase(vector<vector<uchar>> &des);
	bool InitKCenter();
	bool Cluster();
	bool UpdateGroupCenter();
	double DistBetweenDes(vector<uchar> &p1, vector<uchar> &p2);
	double RealDistBetweenDes(vector<uchar> &p1, vector<uchar> &p2);
	bool ExistCenterShift(vector<vector<uchar>> &prev_center,vector<vector<uchar>> &cur_center);
	bool SaveFile(FILE*, FILE*);

	//vector<vector<vector<uchar>>> m_center_seq;

};
#endif // !BASE_H
