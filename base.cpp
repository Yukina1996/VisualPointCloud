#include "base.h"



CLSEstimator::CLSEstimator()
{
	m_nState = 0;
	m_nMeas = 0;
	m_bPDiagonal = false;
	m_bInit = false;
	m_Sigma = 1.0;
	m_SigmaPre = 1.0;
}




CLSEstimator::~CLSEstimator()
{
	FreeLSQ();
}


void CLSEstimator::FreeLSQ()
{
	m_StateX.freeMatrix();
	m_StateP.freeMatrix();
	m_B.freeMatrix();
	m_P.freeMatrix();
	m_L.freeMatrix();
	m_W.freeMatrix();
	m_V.freeMatrix();
}


void CLSEstimator::ClearLSQ()
{
	m_B.zero();
	m_L.zero();
	m_P.zero();
}


bool CLSEstimator::initLSQ(int nState, int nMeas, bool bPDiagonal)
{
	if (nState == m_nState && nMeas == m_nMeas)
	{
		ClearLSQ();
	}
	else
	{
		FreeLSQ();

		m_StateX.mat(nState, 1);
		m_StateP.mat(nState, nState);
		m_B.mat(nMeas, nState);
		m_P.mat(nMeas, nMeas);
		m_L.mat(nMeas, 1);
		m_W.mat(nState, 1);
		m_V.mat(nMeas, 1);
	}

	m_nState = nState;
	m_nMeas = nMeas;
	m_bPDiagonal = bPDiagonal;
	m_bInit = true;

	return true;
}


bool CLSEstimator::LSQ(bool bPostEst)
{
	if (m_bInit == false)
	{
		printf("Forget to Init LSQ\n");
		return false;
	}

	MatrixMultiply_BTPB(m_nMeas, m_nState, m_B.x, m_P.x, m_bPDiagonal, m_StateP.x);
	MatrixMultiply_BTPL(m_nMeas, m_nState, m_B.x, m_P.x, m_L.x, m_bPDiagonal, m_W.x);

	if (m_StateP.inv_sp() == false)
	{
		printf("fail to perform LS, the Co-factor matrix is singular\n");
		return false;
	}

	m_StateX.multiply(m_StateP, m_W);

	m_V.multiply(m_B, m_StateX);
	m_V.subtract_r(m_L); // V = B*X - L

	if (bPostEst)
	{
		m_Sigma = 1.0;

		if (m_B.row > m_B.col)
		{
			MatrixMultiply_BTPB(m_V.row, m_V.col, m_V.x, m_P.x, m_bPDiagonal, &m_Sigma);
			m_Sigma = m_Sigma / (m_B.row - m_B.col);
			if (m_Sigma < 1.0) m_Sigma = 1.0; // sigma不要小于1.0
			m_StateP.scale(m_Sigma);
			m_P.scale(1.0 / m_Sigma);
		}
		else
		{
			m_StateP.scale(m_SigmaPre);
			m_P.scale(1.0 / m_SigmaPre);
		}

		if (m_B.row > m_B.col) m_SigmaPre = m_Sigma;
		m_Sigma = xsqrt(m_SigmaPre);
	}

	return true;
}



CLSEstGeneral::CLSEstGeneral()
{
	m_nStateT = 0; 
	m_DimsT = 0; 
	m_nStateX = 0;
	m_DimsX = 0; 
	m_nMeas = 0; 
	m_DimsMeas = 0;
}

CLSEstGeneral::~CLSEstGeneral()
{
	FreeLSQ();
}

bool CLSEstGeneral::initLSQ(int nStateT, int DimsT, int nStateX, int DimsX, int nMeas, int DimsMeas)
{
	if (nStateT == m_nStateT && DimsT == m_DimsT && nStateX == m_nStateX && DimsX == m_DimsX && nMeas == m_nMeas && DimsMeas == m_DimsMeas)
	{
		ClearLSQ();
	}
	else
	{
		FreeLSQ();
		m_N11.mat(nStateT*DimsT, nStateT*DimsT);
		m_N11i.mat(nStateT*DimsT, nStateT*DimsT);
		m_N22.mat(nStateX*DimsX, nStateX*DimsX);
		m_N12.mat(nStateT*DimsT, nStateX*DimsX);
		m_W1.mat(nStateT*DimsT, 1);
		m_W2.mat(nStateX*DimsX, 1);
		m_N12TN11iN12.mat(nStateX*DimsX, nStateX*DimsX);
		m_N12TN11iW1.mat(nStateX*DimsX, 1);

		m_StateT.mat(nStateT*DimsT, 1);
		m_StateTP.mat(nStateT*DimsT, nStateT*DimsT);
		m_StateX.mat(nStateX*DimsX, 1);
		m_StateXP.mat(nStateX*DimsX, nStateX*DimsX);

	}

	m_nStateT = nStateT;
	m_DimsT = DimsT;
	m_nStateX = nStateX;
	m_DimsX = DimsX;
	m_nMeas = nMeas;
	m_DimsMeas = DimsMeas;

	return true;
}


//和静态变量对应的N22 没有涉及到求逆操作 所以先验信息的权可以最后加
//和动态变量对应的N11 由于每一块N11i都有对应的求逆操作，所以先验信息的权不能再最后加
bool CLSEstGeneral::AddInfo(int id, Mat A, Mat B, Mat L, Mat P, Mat PriP_T, Mat PriLT)
{
	Mat N11i(m_DimsT, m_DimsT), AT(m_DimsT, m_DimsMeas), ATP(m_DimsT, m_DimsMeas);
	Mat N22i(m_nStateX*m_DimsX, m_nStateX*m_DimsX), BT(m_nStateX*m_DimsX, m_DimsMeas), BTP(m_nStateX*m_DimsX, m_DimsMeas);
	Mat N12i(m_DimsT, m_nStateX*m_DimsX), W1i(m_DimsT, 1), W2i(m_nStateX*m_DimsX, 1),N11i_i(m_DimsT, m_DimsT);
	Mat N12Ti(m_nStateX*m_DimsX, m_DimsT),tmp(m_nStateX*m_DimsX, m_DimsT),N12TN11IN12_i(m_nStateX*m_DimsX, m_nStateX*m_DimsX), N12TN11iW1_i(m_nStateX*m_DimsX, 1);
	Mat tmpPri(m_DimsT, 1);
	AT.transpose(A);
	ATP.multiply(AT, P);
	N11i.multiply(ATP, A);
	N11i.add2(PriP_T);
	BT.transpose(B);
	BTP.multiply(BT, P);
	N22i.multiply(BTP, B);
	N12i.multiply(ATP, B);

	W1i.multiply(ATP, L);	
	tmpPri.multiply(PriP_T, PriLT);
	W1i.add2(tmpPri);
	W2i.multiply(BTP, L);

	m_N22.add2(N22i);
	m_W2.add2(W2i);

	N11i_i = N11i; N11i_i.inv();
	N12Ti.transpose(N12i);
	m_N11.setSub(id*m_DimsT, id*m_DimsT, m_DimsT, m_DimsT, N11i);
	m_N11i.setSub(id*m_DimsT, id*m_DimsT, m_DimsT, m_DimsT, N11i_i);
	m_W1.setSub(id*m_DimsT, 0, m_DimsT, 1, W1i);
	m_N12.setSub(id*m_DimsT, 0, m_DimsT, m_nStateX*m_DimsX, N12i);

	tmp.multiply(N12Ti, N11i_i);

	N12TN11IN12_i.multiply(tmp, N12i);
	N12TN11iW1_i.multiply(tmp, W1i);

	m_N12TN11iN12.add2(N12TN11IN12_i);
	m_N12TN11iW1.add2(N12TN11iW1_i);

	
	N11i.freeMatrix();
	AT.freeMatrix();
	ATP.freeMatrix();
	N22i.freeMatrix();
	BT.freeMatrix();
	BTP.freeMatrix();
	N12i.freeMatrix();
	W1i.freeMatrix();
	W2i.freeMatrix();
	N11i_i.freeMatrix();
	N12Ti.freeMatrix();
	tmp.freeMatrix();
	N12TN11IN12_i.freeMatrix();
	N12TN11iW1_i.freeMatrix();
	tmpPri.freeMatrix();

	return true;
}


bool CLSEstGeneral::AddPriInfo(Mat PriP_X, Mat PriLx)
{
	if ((PriP_X.row != m_N22.row) || (PriP_X.col != m_N22.col) || (PriLx.row != m_W2.row))
		return false;
	

	Mat tmp(m_nStateX*m_DimsX, 1);
	m_N22.add2(PriP_X);
	tmp.multiply(PriP_X, PriLx);
	m_W2.add2(tmp);

	tmp.freeMatrix();

	return true;
}

bool CLSEstGeneral::LSQ()
{
	Mat tmpInnoX(m_nStateX*m_DimsX, 1);
	//1.求解静态变量及其方差
	MatrixSubtraction(m_nStateX*m_DimsX, m_nStateX*m_DimsX, m_N22.x, m_N12TN11iN12.x, m_StateXP.x);
	m_StateXP.inv();
	MatrixSubtraction(m_nStateX*m_DimsX, 1, m_W2.x, m_N12TN11iW1.x, tmpInnoX.x);
	MatrixMultiply(m_nStateX*m_DimsX, m_nStateX*m_DimsX, m_StateXP.x, m_nStateX*m_DimsX, 1, tmpInnoX.x, m_StateX.x);
	//2.求解动态变量及其协方差
	for (int i = 0; i < m_nStateT; i++)
	{
		Mat N11i_i(m_DimsT, m_DimsT), W1_i(m_DimsT, 1), N12_i(m_DimsT, m_nStateX*m_DimsX), N12T_i(m_nStateX*m_DimsX, m_DimsT), tmp(m_DimsT, 1);
		Mat tmpInno(m_DimsT, 1), Ti(m_DimsT, 1);
		//2.1 求解动态变量
		m_N11i.getSub(i*m_DimsT, i*m_DimsT, m_DimsT, m_DimsT, N11i_i);
		m_W1.getSub(i*m_DimsT, 0, m_DimsT, 1,W1_i);
		m_N12.getSub(i*m_DimsT, 0, m_DimsT, m_nStateX*m_DimsX, N12_i);
		N12T_i.transpose(N12_i);
		MatrixMultiply(m_DimsT, m_nStateX*m_DimsX, N12_i.x, m_nStateX*m_DimsX, 1, m_StateX.x, tmp.x);
		MatrixSubtraction(m_DimsT, 1, W1_i.x, tmp.x, tmpInno.x);
		MatrixMultiply(m_DimsT, m_DimsT, N11i_i.x, m_DimsT, 1, tmpInno.x, Ti.x);
		m_StateT.setSub(i*m_DimsT, 0, m_DimsT, 1, Ti);
		//2.2 求解动态变量协方差
		Mat tmp1(m_DimsT, m_nStateX*m_DimsX), tmp2(m_DimsT, m_nStateX*m_DimsX), tmp3(m_DimsT, m_DimsT), TPi(m_DimsT, m_DimsT);
		tmp1.multiply(N11i_i, N12_i);
		tmp2.multiply(tmp1, m_StateXP);
		tmp3.multiply(tmp2, N12T_i);
		TPi.multiply(tmp3, N11i_i);
		TPi.add2(N11i_i);
		m_StateTP.setSub(i*m_DimsT, i*m_DimsT, m_DimsT, m_DimsT, TPi);

		N11i_i.freeMatrix();
		W1_i.freeMatrix();
		N12_i.freeMatrix();
		N12T_i.freeMatrix();
		tmp.freeMatrix();
		tmpInno.freeMatrix();
		Ti.freeMatrix();
		tmp1.freeMatrix();
		tmp2.freeMatrix();
		tmp3.freeMatrix();
		TPi.freeMatrix();
	}
	tmpInnoX.freeMatrix();
	return true;
}

void CLSEstGeneral::FreeLSQ()
{
	m_N11.freeMatrix();
	m_N11i.freeMatrix();
	m_N12.freeMatrix();
	m_N12TN11iN12.freeMatrix();
	m_N12TN11iW1.freeMatrix();
	m_N22.freeMatrix();
	m_W1.freeMatrix();
	m_W2.freeMatrix();
	m_StateT.freeMatrix();
	m_StateTP.freeMatrix();
	m_StateX.freeMatrix();
	m_StateXP.freeMatrix();
}

void CLSEstGeneral::ClearLSQ()
{
	m_N11.zero();
	m_N11i.zero();
	m_N12.zero();
	m_N12TN11iN12.zero();
	m_N12TN11iW1.zero();
	m_N22.zero();
	m_W1.zero();
	m_W2.zero();
	m_StateT.zero();
	m_StateTP.zero();
	m_StateX.zero();
	m_StateXP.zero();
}

CKMeans::CKMeans()
{
}

CKMeans::~CKMeans()
{

}

//输入待聚类的特征子数据
bool CKMeans::SetInputDesDatabase(vector<vector<uchar>>& des)
{
	int num_des = des.size();
	m_desDB.resize(num_des);
	for (int i = 0; i < num_des; i++)
	{
		m_desDB[i].des = des[i];
		
		m_desDB[i].GroupID = 0;

	}
	return true;
}

//初始化最初的K个类的中心
bool CKMeans::InitKCenter()
{
	if (m_k == 0)
	{
		printf("在此之前必须要调用setK()函数\n");
		return false;
	}
	
	m_center.resize(m_k);
	for (int i = 0; i < m_k; i++)
	{
		//???????????初值
		int n = m_desDB.size();
		int index = (n / m_k)*i;
	
		for (int j = 0; j < 128; j++)
		{
			m_center[i].push_back(m_desDB[index].des[j]);
		}
		
	}

	return true;
}

//聚类
bool CKMeans::Cluster()
{
	int times = 0;			//计数，以防聚类失败，迭代无穷
	vector<vector<uchar>> v_center(m_center.size());

	do 
	{
		for (int i = 0; i < m_desDB.size(); i++)
		{
			double min_dist = DBL_MAX;
			int grp_id = 0;
			for (int j = 0; j < m_k; j++)
			{
				double dist = DistBetweenDes(m_desDB[i].des, m_center[j]);
				if (min_dist - dist > 0.000001)
				{
					min_dist = dist;
					grp_id = j;
				}
			}
			m_grp_desDB[grp_id].push_back(Des_t(m_desDB[i].des, grp_id));
		}


		//保存上一次迭代的中心店
		for (int i = 0; i < m_center.size(); i++)
		{
			v_center[i] = m_center[i];
		}

		if (!UpdateGroupCenter())
		{
			return false;
		}


		if (!ExistCenterShift(v_center, m_center))
		{
			break;
		}

		if (times >= 50)
		{
			break;
		}

		for (int i = 0; i < m_grp_desDB.size(); i++)
		{
			m_grp_desDB[i].clear();
		}

		times++;
	} while (true);



	return true;
}

//更新K类的中心??
bool CKMeans::UpdateGroupCenter()
{
	if (m_center.size() != m_k)
	{
		printf("类别的个数不为K\n");
		return false;
	}

	for (int i = 0; i < m_k; i++)
	{
		int des_num_in_grp = m_grp_desDB[i].size();
		vector<uchar> des_mean(0);
		des_mean.resize(128);
		for (int j = 0; j < des_num_in_grp; j++)
		{
			for (int u = 0;  u < 128;  u++)
			{
				des_mean[u] = des_mean[u] + m_grp_desDB[i][j].des[u] / des_num_in_grp;
			}
		}
		m_center[i] = des_mean;
	}
	
	return true;
}

//计算两个特征子间的欧氏距离
double CKMeans::DistBetweenDes(vector<uchar>& p1, vector<uchar>& p2)
{
	double dist = 0.0;
	for (int i = 0; i < 128; i++)
	{
		dist = dist + (p1[i] - p2[i])*(p1[i] - p2[i]);
	}
	dist = sqrt(dist);

	return dist;
}

//计算两个特征子间真实的的欧氏距离，去归一化
double CKMeans::RealDistBetweenDes(vector<uchar>& p1, vector<uchar>& p2)
{
	double dist = 0.0;
	double a = 0.0;
	double b = 0.0;
	for (int i = 0; i < 128; i++)
	{
		a = p1[i] / 512.0;
		b = p2[i] / 512.0;
		dist = dist + (a - b)*(a - b);
	}
	dist = sqrt(dist);

	return dist;
}

//是否存在中心点移动
bool CKMeans::ExistCenterShift(vector<vector<uchar>>& prev_center, vector<vector<uchar>>& cur_center)
{
	for (int i = 0; i < m_k; i++)
	{
		double dist = DistBetweenDes(prev_center[i], cur_center[i]);
		if (dist > 0.001)
		{
			return true;
		}
	}

	return false;
}

//将聚类的点分别存到des的bin文件和txt文件中??
bool CKMeans::SaveFile(FILE* fdes, FILE* fdes2)
{
	//for (int i = 0; i < m_k; i++)
	//{
	//	fwrite(&m_center[i], sizeof(uchar), 128, fdes);
	//	for (int j = 0; j < 128; j++)
	//	{
	//		fprintf(fdes2, "%4u", m_center[i][j]);
	//	}
	//	fprintf(fdes2, "\n");
	//}
	
	return true;
}



