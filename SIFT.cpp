#include "SIFT.h"

CPointCloud::CPointCloud()
{
	m_fx = 0.0;
	m_fy = 0.0;
	m_cx = 0.0;
	m_cy = 0.0;
	m_bf = 0.0;
	m_VarPix = 0.0;
	m_BackProjectUVError = 0.0;
	m_BAState = 0;
	m_PicFolderPath = "";
	m_DesFilePath = "";
	m_PointFilePath = "";
	m_UVFilePath = "";
	m_PointInfoPath = "";
	m_PointFileTxt = "";
	m_DesFileTxt = "";


	m_frameIndex = 0;
	m_DesIndex = 0;
	m_PointIndex = 0;
	m_ValidDesIndex = 0;
	m_ValidPointIndex = 0;

	m_CurFrameAtt[9] = { 0.0 };  
	m_CurFramePos[3] = { 0.0 };
}

CPointCloud::~CPointCloud()
{

}

bool CPointCloud::runPointCloud(string optfile, CPCOPTION * pPCOpt)
{
	if (pPCOpt) m_PCOpt = *pPCOpt;
	if (m_PCOpt.readOptFile(optfile) == false) return false;
	//1.初始化
	if (Init(&m_PCOpt) == false) return false;

	//2.从序列影像和真值位姿中增量构建特征点的三维地图数据
	for (int j = 0; j < m_vTimestamps.size(); j++)
	{
		m_frameIndex = j;


		//获取所有影像的左右图像匹配情况
		if (RunOneLRSIFT() == false)
			return false;

		RunSeqMatch();
		InitPoint();
		UpdateInfo();
	}

	//3.将构建出的每个地图点采用不同的像素观测值进行平差
 	DealAllPoint();	

	//4.筛选地图数据库中的地图点并确定其对应的特征描述子们
    //SelectPoint();
    //5.保存剩余地图点文件
    //SaveFile();

	return true;
}


bool CPointCloud::Init(CPCOPTION* pPCOpt)
{
	//1.基本参数赋值
	m_fx = pPCOpt->m_fx;
	m_fy = pPCOpt->m_fy;
	m_cx = pPCOpt->m_cx;
	m_cy = pPCOpt->m_cy;
	m_bf = pPCOpt->m_bf;
	m_PicFolderPath = pPCOpt->m_PICFolder;
	m_DesFilePath = pPCOpt->m_DesFilePath;
	m_PointFilePath = pPCOpt->m_PointFilePath;
	m_UVFilePath = pPCOpt->m_UVFilePath;
	m_PointInfoPath = pPCOpt->m_PointInfoPath;
	m_PointFileTxt = m_PicFolderPath + "\\Point.txt";
	m_DesFileTxt = m_PicFolderPath + "\\Des.txt";
	string SiftFolder = m_PicFolderPath + "\\SIFT";
	_mkdir(SiftFolder.c_str());

	m_VarPix = pPCOpt->m_VarPix;
	m_BackProjectUVError = pPCOpt->m_BackProjectUVError;
	//2.时间和图片路径
	ifstream fTimes;
	fTimes.open(pPCOpt->m_TimesPath.c_str());
	double t = 0;
	int week = 0;
	GPSTIME gt;
	while (!fTimes.eof())
	{
		string s;
		getline(fTimes, s);
		if (!s.empty())
		{
			stringstream ss;
			ss << s;
			ss >> t >> week;

			if (week == 0)
			{
				week = 999;
			}
			gt.GPSWeek = week;
			gt.secsOfWeek = int(t);
			gt.fracOfSec = t - gt.secsOfWeek;
			m_vTimestamps.push_back(gt);
		}
	}
	const int nTimes = m_vTimestamps.size();
	m_vstrImageLeft.resize(nTimes);
	m_vstrImageRight.resize(nTimes);

	//提取后缀名
	string Suffix;
	if (m_PCOpt.m_SuffixFlag == 0)
		Suffix = ".jpg";
	else if (m_PCOpt.m_SuffixFlag == 1)
		Suffix = ".png";
	else
		Suffix = ".bmp";

	for (int i = 0; i < nTimes; i++)
	{
		stringstream ss;
		ss << setfill('0') << setw(6) << i + m_PCOpt.m_StartFid;
		m_vstrImageLeft[i] = pPCOpt->m_PicLFolderPath + ss.str() + Suffix;
		m_vstrImageRight[i] = pPCOpt->m_PicRFolderPath + ss.str() + Suffix;
	}

	//3.真值位置、姿态读取
	if(readGTFile(pPCOpt->m_GTFilePath)==false)
		return false;

	return true;
}

bool CPointCloud::readGTFile(string path)
{
	FILE* fgt = fopen(path.c_str(), "rt");
	if (fgt == NULL) {
		printf("Error:unable to open groundtruth file!\n");
		return false;
	}

	int num = 0;
	while (!feof(fgt))
	{
		GroundTruth pose;
		double tmp;
		fscanf(fgt, "%d %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &pose.gt.GPSWeek, &tmp, &pose.AmbState, &pose.Pos[0], &pose.Pos[1], &pose.Pos[2], &pose.VarPos[0], &pose.VarPos[4], &pose.VarPos[8], &pose.VarPos[1], &pose.VarPos[2], &pose.VarPos[5], &pose.Att[0], &pose.Att[1], &pose.Att[2], &pose.VarAtt[0], &pose.VarAtt[4], &pose.VarAtt[8], &pose.VarAtt[1], &pose.VarAtt[2], &pose.VarAtt[5]);
		pose.gt.secsOfWeek = int(tmp);
		pose.gt.fracOfSec = tmp - pose.gt.secsOfWeek;
		RotaAngle2RotaMatrix(pose.Att, pose.Rwc);
		RotaMatrix2IMatrix(pose.Rwc);
		MatrixCopy(3, 3, pose.VarAtt, pose.PAtt);
		MatrixInv(3, 3, pose.PAtt);
		MatrixCopy(3, 3, pose.VarPos, pose.PPos);
		MatrixInv(3, 3, pose.PPos);
		m_vGroundTruth.push_back(pose);
		num++;
	}
	fclose(fgt);

	if (num != m_vTimestamps.size()) return false;

	return true;
}

bool CPointCloud::RunOneLRSIFT()
{
	string ImgLeftPath = m_vstrImageLeft[m_frameIndex];
	string ImgRightPath = m_vstrImageRight[m_frameIndex];

	printf("%s\n", ImgLeftPath.c_str());
	SiftGPU sift;
	SiftMatchGPU matcher;
	char* argv[] = { (char*)"-fo",(char*)"-1",(char*)"-v" ,(char*)"0"};
	int argc = sizeof(argv) / sizeof(char*);
	sift.ParseParam(argc, argv);
	if (sift.CreateContextGL() != SiftGPU::SIFTGPU_FULL_SUPPORTED)
		return false;

	int numl, numr;
	vector<SiftGPU::SiftKeypoint> keysl(1), keysr(1);
	vector<float> desl(1), desr(1);

	//2.左右影像特征提取
	sift.RunSIFT(ImgLeftPath.c_str());
	numl = sift.GetFeatureNum();
	keysl.resize(numl);
	desl.resize(numl * 128);
	sift.GetFeatureVector(&keysl[0], &desl[0]);

	sift.RunSIFT(ImgRightPath.c_str());
	numr = sift.GetFeatureNum();
	keysr.resize(numr);
	desr.resize(numr * 128);
	sift.GetFeatureVector(&keysr[0], &desr[0]);
	
	//3.左右影像匹配

	matcher.SetLanguage(SiftMatchGPU::SIFTMATCH_GLSL);
	//matcher.SetLanguage(SiftMatchGPU::SIFTMATCH_CUDA);
	matcher.VerifyContextGL();
	matcher.SetDescriptors(0, numl, &desl[0]);
	matcher.SetDescriptors(1, numr, &desr[0]);
	int(*match_buf)[2] = new int[numl][2];
	int num_match = matcher.GetSiftMatch(numl, match_buf);

	//FILE* fp = fopen("F:\\KITTIData\\dataset\\sequences\\00\\depth_wyy.txt", "a+");
	//float z_depth;

	//4.提取匹配的点对和描述子
	for (int k = 0; k < num_match; ++k)
	{
		SiftGPU::SiftKeypoint & key1 = keysl[match_buf[k][0]];
		SiftGPU::SiftKeypoint & key2 = keysr[match_buf[k][1]];
		//z_depth = m_bf / (key1.x - key2.x);
		//fprintf(fp, "%f\n", z_depth);

		if ((abs(key1.y - key2.y) > 1) || (key1.x < key2.x) || (abs(m_bf / (key1.x - key2.x)) > (100 * m_bf / m_fx)))
			continue;
		
		//if ((key1.y > 1040) || (key2.y > 1040))
		//	continue;
		//if ((key1.y > 1075) || (key2.y > 1075))
		//	continue;
		
		m_Cursiftinfo.Keys.push_back(key1);
		m_Cursiftinfo.UR.push_back(key2.x);
	
		for (int i = 0; i < 128; i++)
		{
			//如果GlobalUtil::_NormalizedSIFT 经过归一化 ，则SIFTGPU->savesift函数存储描述符时，会做如下还原处理
			//floor(0.5+512.0f*(*pd))
			//而在这存储的是归一化后的值
			unsigned int a = ((unsigned int)floor(0.5 + 512.0f*(desl[128 * match_buf[k][0] + i])));
			m_Cursiftinfo.des.push_back(uchar(a));
		}
	}
	m_Cursiftinfo.KeyNum = m_Cursiftinfo.Keys.size();
	m_Cursiftinfo.Point_index.resize(m_Cursiftinfo.KeyNum);

	//fclose(fp);
	free(match_buf);
	return true;
}


void CPointCloud::RunSeqMatch()
{
	//1.SIFTGPU参数设置
	SiftMatchGPU matcher;
	//matcher.SetLanguage(SiftMatchGPU::SIFTMATCH_CUDA);
	matcher.SetLanguage(SiftMatchGPU::SIFTMATCH_GLSL);
	matcher.VerifyContextGL();

	//FILE* fp = fopen("F:\\KITTIData\\dataset\\sequences\\00\\epipolar_wyy.txt", "a+");


	if (m_frameIndex < 20)
	{
		for (int i = m_frameIndex - 1; i > -1; i--)
		{
			//3.左右影像匹配
	        //代码中并未采用欧式距离
	        //而是将128维的特征向量归一化，点乘两个特征向量求acos，阈值设为0.8,最优和次优比值默认为0.7
			SIFTInfo sinfo;
			ReadSiftInfo(i, &sinfo);
			matcher.SetDescriptors(0, sinfo.KeyNum, &sinfo.des[0]);
			matcher.SetDescriptors(1, m_Cursiftinfo.KeyNum, &m_Cursiftinfo.des[0]);
			int(*match_buf)[2] = new int[sinfo.KeyNum][2];
			int num_match = matcher.GetSiftMatch(sinfo.KeyNum, match_buf);

			for (int k = 0; k < num_match; ++k)
			{
				int LastIndex = match_buf[k][0];
				int CurIndex = match_buf[k][1];
				if (m_Cursiftinfo.Point_index[CurIndex] == 0)
				{
					double Cur_uv[3] = { m_Cursiftinfo.Keys[CurIndex].x,m_Cursiftinfo.Keys[CurIndex].y,1.0 };
					double Match_uv[3] = { sinfo.Keys[LastIndex].x,sinfo.Keys[LastIndex].y,1.0 };
					double result = EssentialVerify(i, Cur_uv, Match_uv);
					if (result < 10)
					{
						m_Cursiftinfo.Point_index[CurIndex] = sinfo.Point_index[LastIndex];
					}
					//fprintf(fp, "%lf\n", result);
					
				}
			}
			free(match_buf);
		}
	}		
	else
	{
		for (int i = m_frameIndex - 1; i > m_frameIndex-21; i--)
		{
			SIFTInfo sinfo;
			ReadSiftInfo(i, &sinfo);

			matcher.SetDescriptors(0, sinfo.KeyNum, &sinfo.des[0]);
			matcher.SetDescriptors(1, m_Cursiftinfo.KeyNum, &m_Cursiftinfo.des[0]);
			int(*match_buf)[2] = new int[sinfo.KeyNum][2];
			int num_match = matcher.GetSiftMatch(sinfo.KeyNum, match_buf);

			//4.提取匹配的点对和描述子
			for (int k = 0; k < num_match; ++k)
			{
				int LastIndex = match_buf[k][0];
				int CurIndex = match_buf[k][1];
				if (m_Cursiftinfo.Point_index[CurIndex] == 0)
				{
					double Cur_uv[3] = { m_Cursiftinfo.Keys[CurIndex].x,m_Cursiftinfo.Keys[CurIndex].y,1.0 };
					double Match_uv[3] = { sinfo.Keys[LastIndex].x,sinfo.Keys[LastIndex].y,1.0 };
					double result = EssentialVerify(i, Cur_uv, Match_uv);
					if (result < 10)
					{
						m_Cursiftinfo.Point_index[CurIndex] = sinfo.Point_index[LastIndex];
					}
					//fprintf(fp, "%lf\n", result);
				}
			}
			free(match_buf);
		}
	}

	//fclose(fp);
}

double CPointCloud::EssentialVerify(int Frameidx, double Cur_uv[], double Match_uv[])
{
	GroundTruth Cur_GT = m_vGroundTruth[m_frameIndex];
	GroundTruth Match_GT = m_vGroundTruth[Frameidx];

	double K[9], iK[9], iKT[9], iKTE[9];
	double result;
	double tmp[3], dP[3];
	double Rc1w[9], R12[9], t[3], tx[9], E[9], F[9];
	//内参矩阵
	K[0] = m_fx; K[1] = 0.0 ; K[2] = m_cx;
	K[3] = 0.0 ; K[4] = m_fy; K[5] = m_cy;
	K[6] = 0.0 ; K[7] = 0.0 ; K[8] = 1   ;
	MatrixCopy(3, 3, K, iK);
	MatrixInv(3, 3, iK);
	MatrixCopy(3, 3, iK, iKT);
	MatrixTranspose(3, 3, iKT);

	//求得当前帧对应的Rcw和T
	MatrixCopy(3, 3, Cur_GT.Rwc, Rc1w);
	MatrixTranspose(3, 3, Rc1w);

	//求本质矩阵R12,t
	MatrixMultiply(3, 3, Match_GT.Rwc, 3, 3, Rc1w, R12);
	MatrixSubtraction(3, 1, Cur_GT.Pos, Match_GT.Pos, dP);
	MatrixMultiply(3, 3, Match_GT.Rwc, 3, 1, dP, t);
	MatrixSkewSymmetric(t, tx);
	MatrixMultiply(3, 3, tx, 3, 3, R12, E);
	MatrixMultiply(3, 3, iKT, 3, 3, E, iKTE);
	MatrixMultiply(3, 3, iKTE, 3, 3, iK, F);

	//MatrixMultiply(1, 3, Cur_uv, 3, 3, F, tmp);
	MatrixMultiply(3, 3, F, 3, 1, Match_uv, tmp);
	
	//result = tmp[0] * Match_uv[0] + tmp[1] * Match_uv[1] + tmp[2] * Match_uv[2];
	result = Cur_uv[0] * tmp[0] + Cur_uv[1] * tmp[1] + Cur_uv[2] * tmp[2];

	result = result / sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);

	return abs(result);
}



bool CPointCloud::InitPoint()
{
	stringstream ss;
	ss << setfill('0') << setw(6) << m_frameIndex;
	string path = m_PicFolderPath + "\\SIFT\\"+ss.str() + ".sift";
	FILE* fdes = fopen(path.c_str(), "wb");
	if (fdes == NULL) {
		printf("Error:unable to open descriptor file!\n");
		return false;
	}
	fwrite(&m_Cursiftinfo.KeyNum, sizeof(int), 1, fdes);
	double Rcw[9];
	for (int j = 0; j < 3; j++)
	{
		for (int k = 0; k < 3; k++)
		{
			m_CurFrameAtt[3 * j + k] = m_vGroundTruth[m_frameIndex].Rwc[3 * j + k];
			Rcw[3 * j + k] = m_vGroundTruth[m_frameIndex].Rwc[3 * j + k];//此时得到的是Rwc；
		}
		m_CurFramePos[j] = m_vGroundTruth[m_frameIndex].Pos[j];
	}

	for (int i = 0; i < m_Cursiftinfo.KeyNum; i++)
	{
		UV uv;
		double xyz[3],XYZ[3];

		uv.Frame_id = m_frameIndex;
		uv.Des_id = m_DesIndex;
		uv.UV[0] = m_Cursiftinfo.Keys[i].x;
		uv.UV[1] = m_Cursiftinfo.Keys[i].y;
		uv.UV[2] = uv.UV[0] - m_Cursiftinfo.UR[i];
		uv.ValidFlag = 1;
		if (m_Cursiftinfo.Point_index[i] == 0)
		{
			m_Cursiftinfo.Point_index[i] = m_PointIndex;
			m_AllPoint[m_PointIndex].uvlist.push_back(uv);
			BacProject(uv.UV, xyz);
			MatrixInv(3, 3, Rcw);
			MatrixMultiply(3, 3, Rcw, 3, 1, xyz, XYZ);
			MatrixAddition2(3, 1, m_CurFramePos, XYZ);
			
			m_AllPoint[m_PointIndex].Point_id = m_PointIndex;
			m_AllPoint[m_PointIndex].n_Frame = 1;
			m_AllPoint[m_PointIndex].XYZ[0] = XYZ[0];
			m_AllPoint[m_PointIndex].XYZ[1] = XYZ[1];
			m_AllPoint[m_PointIndex].XYZ[2] = XYZ[2];

			m_PointIndex++;
		}
		else
		{
			int p_idx = m_Cursiftinfo.Point_index[i];
			int Repeat_flag = 0;
			for (int k1 = 0; k1 < m_AllPoint[p_idx].uvlist.size(); k1++)
			{
				if (m_AllPoint[p_idx].uvlist[k1].Frame_id == uv.Frame_id)
				{
					Repeat_flag = 1;
					break;
				}
			}
			if (Repeat_flag == 0)
			{
				m_AllPoint[p_idx].uvlist.push_back(uv);
				m_AllPoint[p_idx].n_Frame++;
			}
		}
		fwrite(&m_DesIndex, sizeof(int), 1, fdes);
		fwrite(&m_Cursiftinfo.Point_index[i], sizeof(int), 1, fdes);
		fwrite(&m_Cursiftinfo.Keys[i].x, sizeof(float), 1, fdes);//为了后续极限约束用到 所以需要把坐标也输出
		fwrite(&m_Cursiftinfo.Keys[i].y, sizeof(float), 1, fdes);
		fwrite(&m_Cursiftinfo.des[i * 128], sizeof(uchar), 128, fdes);
		m_DesIndex++;
	}

	fclose(fdes);
}



void CPointCloud::UpdateInfo()
{
	m_Cursiftinfo.des.resize(0);
	m_Cursiftinfo.Keys.resize(0);
	m_Cursiftinfo.UR.resize(0);
	m_Cursiftinfo.Point_index.resize(0);
	m_Cursiftinfo.KeyNum = 0;
}


int CPointCloud::BAPoint(int Point_id)
{
	CurPoint = m_AllPoint[Point_id];
	while (CurPoint.n_Frame>1)
	{
		//1.做一次BA
		//每一次BA开始前，状态置为1；
		m_BAState = 1;
		RunOneBA();
		//2.计算像素的反投影误差 判断反投影误差是否超过阈值 
		CalBackProjectError();
		if (m_BAState == 1)
		{
			printf("%d th point:BA success\n",Point_id);
			break;
		}
		else
		{
			//3.查过残差设定阈值，找到残差最大的对应观测值，并剔除
			int MeasID = 0;
			double MaxMeanBacProjectError = 0.0;
			for (int j = 0; j < CurPoint.uvlist.size(); j++)
			{
				if (CurPoint.uvlist[j].ValidFlag > 0 && CurPoint.uvlist[j].MeanBP_UVError > MaxMeanBacProjectError)
				{
					MeasID = j;
				}
			}
			CurPoint.uvlist[MeasID].ValidFlag = 0;
			CurPoint.n_Frame = CurPoint.n_Frame - 1;
			//CurPoint的初值变了，保持初值不变
			CurPoint.XYZ[0] = m_AllPoint[Point_id].XYZ[0];
			CurPoint.XYZ[1] = m_AllPoint[Point_id].XYZ[1];
			CurPoint.XYZ[2] = m_AllPoint[Point_id].XYZ[2];
		}
	}

	//BA_State最终为1 表示平差成功
	if (m_BAState == 1)
	{
		//进一步判断
		//这个判断还有待考虑
		if (CurPoint.VarXYZ[0] < 4 && CurPoint.VarXYZ[4] < 4 && CurPoint.VarXYZ[8] < 4)
		{
			//更新m_AllPoint的信息
			m_AllPoint[Point_id].XYZ[0] = CurPoint.XYZ[0];
			m_AllPoint[Point_id].XYZ[1] = CurPoint.XYZ[1];
			m_AllPoint[Point_id].XYZ[2] = CurPoint.XYZ[2];
			for (int i = 0; i < 9; i++)
				m_AllPoint[Point_id].VarXYZ[i] = CurPoint.VarXYZ[i];
			m_AllPoint[Point_id].n_Frame = CurPoint.n_Frame;
			m_AllPoint[Point_id].uvlist.resize(0);
			int RIF = 0;
			for (int j = 0; j < CurPoint.uvlist.size(); j++)
			{
				if (CurPoint.uvlist[j].ValidFlag > 0)
				{
					m_AllPoint[Point_id].uvlist.push_back(CurPoint.uvlist[j]);
					if (CurPoint.uvlist[j].Frame_id > RIF)
						RIF = CurPoint.uvlist[j].Frame_id;
				}
			}
			m_AllPoint[Point_id].Recent_FI = RIF;
			return 2;
		}
		else
		{
			return 1;
		}
	}
	else
		return 1;
}





int CPointCloud::RunOneBA()
{
	int it = 0;
	double norm1 = 0, norm2 = 0, Di = 0;
	CLSEstGeneral LSEG;

	int validNum = 0;
	for (int i = 0; i < CurPoint.uvlist.size(); i++)
	{
		if (CurPoint.uvlist[i].ValidFlag == 1)
			validNum++;
	}

	

	while (it++ < 10)
	{
		int idx = 0;
		LSEG.initLSQ(validNum, 6, 1, 3, validNum, 3);
		for (int i = 0; i < CurPoint.uvlist.size(); i++)
		{
			if (!CurPoint.uvlist[i].ValidFlag)
				continue;

			int F_id = CurPoint.uvlist[i].Frame_id;
			double P_tmp[3], P_tmpX[9], p[3];
			double du_dp[9], dp_dphi[9], dp_dPc[9], dp_dPp[9];
			double du_dphi[9], du_dPc[9], du_dPp[9], uvd[3], l[3];
			MatrixSubtraction(3, 1, CurPoint.XYZ, m_vGroundTruth[F_id].Pos, P_tmp);
			MatrixSkewSymmetric(P_tmp, P_tmpX);
			MatrixMultiply(3, 3, m_vGroundTruth[F_id].Rwc, 3, 1, P_tmp, p);
			//du_dp
			du_dp[0] = m_fx / p[2]; du_dp[1] = 0; du_dp[2] = -1.0*m_fx*p[0] / (p[2] * p[2]);
			du_dp[3] = 0; du_dp[4] = m_fy / p[2]; du_dp[5] = -1.0*m_fy*p[1] / (p[2] * p[2]);
			du_dp[6] = 0; du_dp[7] = 0; du_dp[8] = -1.0*m_bf / (p[2] * p[2]);
			//dp_dphi
			MatrixMultiply(3, 3, m_vGroundTruth[F_id].Rwc, 3, 3, P_tmpX, dp_dphi);
			M33Scale(-1.0, dp_dphi);
			//dp_dPc
			M33EQU(m_vGroundTruth[F_id].Rwc, dp_dPc);
			M33Scale(-1.0, dp_dPc);
			//dp_dPp
			M33EQU(m_vGroundTruth[F_id].Rwc, dp_dPp);
			//du_dphi,du_dPc,du_dPp
			MatrixMultiply(3, 3, du_dp, 3, 3, dp_dphi, du_dphi);
			MatrixMultiply(3, 3, du_dp, 3, 3, dp_dPc, du_dPc);
			MatrixMultiply(3, 3, du_dp, 3, 3, dp_dPp, du_dPp);
			//残差
			StereoProject(p, uvd);
			MatrixSubtraction(3, 1, CurPoint.uvlist[i].UV, uvd, l);

			//A B  L P 矩阵赋值
			Mat A(3, 6), B(3, 3), L(3, 1), P(3, 3);
			Mat PriP_T(6, 6), PriLT(6, 1);
			A.setSub(0, 0, 3, 3, du_dphi);
			A.setSub(0, 3, 3, 3, du_dPc);
			B.setSub(0, 0, 3, 3, du_dPp);
			L.setSub(0, 0, 3, 1, l);
			P.x[_MI(P, 0, 0)] = 1 / m_VarPix;
			P.x[_MI(P, 1, 1)] = 1 / m_VarPix;
			P.x[_MI(P, 2, 2)] = 1 / m_VarPix;
			PriP_T.setSub(0, 0, 3, 3, m_vGroundTruth[F_id].PAtt);
			PriP_T.setSub(3, 3, 3, 3, m_vGroundTruth[F_id].PPos);
			LSEG.AddInfo(idx, A, B, L, P, PriP_T, PriLT);
			A.freeMatrix();
			B.freeMatrix();
			L.freeMatrix();
			P.freeMatrix();
			PriP_T.freeMatrix();
			PriLT.freeMatrix();
			idx++;
		}
		//解算
		LSEG.LSQ();
		MatrixAddition2(3, 1, LSEG.m_StateX.x, CurPoint.XYZ);
		norm1 = LSEG.m_StateX.norm2();

		if (norm1 < 1.0e-4)
		{
			break;
		}
		if (norm1 - norm2 > 0)
		{
			Di++;
			if (Di > 1)
				break;
		}
		norm2 = norm1;
	}
	

	for (int i = 0; i < 9; i++)
	{
		CurPoint.VarXYZ[i] = LSEG.m_StateXP.x[i];
	}

	return 1;
}

bool CPointCloud::DealAllPoint()
{

	FILE* fout_uv = fopen(m_UVFilePath.c_str(), "w");
	if (fout_uv == NULL) {
		printf("Error:unable to open UVFilePath file!\n");
		return false;
	}

	FILE* fout_pi = fopen(m_PointInfoPath.c_str(), "w");
	if (fout_pi == NULL)
	{
		printf("Error:unable to open PointInfoPath file!\n");
		return false;
	}

	FILE* fp;
	if ((fp = fopen(m_PointFilePath.c_str(), "wb")) == NULL)
	{
		printf("can't open fp file\n");
		return false;
	}

	FILE* fdes;
	if ((fdes = fopen(m_DesFilePath.c_str(), "wb")) == NULL)
	{
		printf("can't open fdes file\n");
		return false;
	}

	FILE* fp2;
	if ((fp2 = fopen(m_PointFileTxt.c_str(), "w")) == NULL)
	{
		printf("can't open fp2 file\n");
		return false;
	}

	FILE* fdes2;
	if ((fdes2 = fopen(m_DesFileTxt.c_str(), "w")) == NULL)
	{
		printf("can't open fdes2 file\n");
		return false;
	}



	fwrite(&m_ValidPointIndex, sizeof(int), 1, fp);
	fwrite(&m_ValidDesIndex, sizeof(int), 1, fp);
	int k1 = 0;
	int k2 = 0;
	for (int idx = 0; idx < m_AllPoint.size(); idx++)
	{
		if (m_AllPoint[idx].n_Frame < 2 || m_AllPoint[idx].n_Frame > 50)
		//if (m_AllPoint[idx].n_Frame < 2 || m_AllPoint[idx].n_Frame > 100)
		{
			continue;
		}
		else
		{
			if (BAPoint(idx) == 2)
			{
				PInfo pinfo;
				PointInfo(fout_uv, fout_pi, idx, &pinfo);
				SelectPointTmp(fp, fdes, fp2, fdes2, idx, &pinfo);
			}
		}


	}
	fseek(fp, 0L, SEEK_SET);
	fseek(fdes, 0L, SEEK_SET);
	fwrite(&m_ValidPointIndex, sizeof(int), 1, fp);
	fwrite(&m_ValidDesIndex, sizeof(int), 1, fp);


	printf("%d %d\n", m_ValidPointIndex, m_ValidDesIndex);

	fclose(fout_uv);
	fclose(fdes);
	fclose(fp);
	fclose(fdes2);
	fclose(fp2);
	return true;
}

void CPointCloud::PointInfo(FILE* fout_uv, FILE* fout_pi, int Point_id, PInfo *pinfo)
{
	vector<double> LateralDis;   //侧向距离
	vector<double> Dis;          //几何距离
	vector<double> Angle;        //视线夹角
	vector<double> Depth;		 //深度距离
	vector<double> Residual;     //u v ul-ur残差和 
	double MMS_Res[3] = { 0.0 }; 
	double HTH[9] = { 0.0 };
	int i;

	PointXYZ P = m_AllPoint[Point_id];
	for (int i = 0; i < P.uvlist.size(); i++)
	{
		double P_tmp[3], P_tmpX[9], p[3], uvd[3], l[3];
		
		int F_id = P.uvlist[i].Frame_id;
		MatrixSubtraction(3, 1, P.XYZ, m_vGroundTruth[F_id].Pos, P_tmp);
		MatrixSkewSymmetric(P_tmp, P_tmpX);
		MatrixMultiply(3, 3, m_vGroundTruth[F_id].Rwc, 3, 1, P_tmp, p);
		StereoProject(p, uvd);
		MatrixSubtraction(3, 1, P.uvlist[i].UV, uvd, l);

		double latDis = sqrt(p[0] * p[0] + p[1] * p[1]);
		double L = sqrt(latDis * latDis + p[2] * p[2]);
		double ang = acos(p[2] / L)*R2D;
		double res = l[0] * l[0] + l[1] * l[1] + l[2] * l[2];
		LateralDis.push_back(latDis);
		Dis.push_back(L);
		Angle.push_back(ang);
		Depth.push_back(abs(p[2]));
		Residual.push_back(res);

		HTH[0] = HTH[0] + m_fx * m_fx / (p[2] * p[2]);
		HTH[2] = HTH[2] - m_fx * m_fx *p[0] / (p[2] * p[2] * p[2]);
		HTH[4] = HTH[4] + m_fy * m_fy / (p[2] * p[2]);
		HTH[5] = HTH[5] - m_fy * m_fy *p[1] / (p[2] * p[2] * p[2]);
		HTH[6] = HTH[2];
		HTH[7] = HTH[5];
		HTH[8] = HTH[8] + (m_fx * m_fx *p[0] * p[0] + m_fy * m_fy *p[1] * p[1] + m_bf * m_bf) / (p[2] * p[2] * p[2] * p[2]);
	}
	GetMedMeanSTD(LateralDis, pinfo->MMS_Lat);
	GetMedMeanSTD(Dis, pinfo->MMS_Dis);
	GetMedMeanSTD(Angle, pinfo->MMS_Ang);
	GetMedMeanSTD(Depth, pinfo->MMS_Dep);
	GetMedMeanSTD(Residual, MMS_Res);
	MatrixInv(3, 3, HTH);
	pinfo->Res = sqrt(MMS_Res[1] / 3);
	pinfo->Post_sigma2 = (MMS_Res[1] * P.uvlist.size() * 1 / m_VarPix) / (3 * P.uvlist.size() - 3);
	pinfo->Dop = MatrixTrace(3, 3, HTH);
	//if (pinfo->MMS_Ang[1] < 90 && pinfo->Res < 10)
	if (pinfo->MMS_Ang[1] < 90)
	{
		fprintf(fout_uv, "%10d %13.4lf %13.4lf %13.4lf ", Point_id, P.XYZ[0], P.XYZ[1], P.XYZ[2]);
		fprintf(fout_uv, "%10d %13.4lf %13.4lf %13.4lf ", P.n_Frame, P.VarXYZ[0] , P.VarXYZ[4] , P.VarXYZ[8]);
		fprintf(fout_uv, "%13.4lf %13.4lf %13.4lf ", pinfo->MMS_Dis[0], pinfo->MMS_Dis[1], pinfo->MMS_Dis[2]);
		fprintf(fout_uv, "%13.4lf %13.4lf %13.4lf ", pinfo->MMS_Lat[0], pinfo->MMS_Lat[1], pinfo->MMS_Lat[2]);
		fprintf(fout_uv, "%13.4lf %13.4lf %13.4lf ", pinfo->MMS_Dep[0], pinfo->MMS_Dep[1], pinfo->MMS_Dep[2]);
		fprintf(fout_uv, "%13.4lf %13.4lf %13.4lf ", pinfo->MMS_Ang[0], pinfo->MMS_Ang[1], pinfo->MMS_Ang[2]);
		fprintf(fout_uv, "%13.4lf %13.4lf %13.4lf\n", pinfo->Res, pinfo->Dop, pinfo->Post_sigma2);
		for (i = 0; i < P.n_Frame; i++)
		{	fprintf(fout_pi, "%13.4lf", Dis[i]);}		fprintf(fout_pi, "\n");
		for (i = 0; i < P.n_Frame; i++)
		{	fprintf(fout_pi, "%13.4lf", LateralDis[i]);}	fprintf(fout_pi, "\n");
		for (i = 0; i < P.n_Frame; i++)
		{	fprintf(fout_pi, "%13.4lf", Depth[i]);}		fprintf(fout_pi, "\n");
		for (i = 0; i < P.n_Frame; i++)
		{	fprintf(fout_pi, "%13.4lf", Angle[i]);}		fprintf(fout_pi, "\n");
		for (i = 0; i < P.n_Frame; i++)
		{	fprintf(fout_pi, "%13.4lf", Residual[i]);}	fprintf(fout_pi, "\n");
	}
}

bool CPointCloud::SelectPoint()
{

	FILE* fpt = fopen(m_PointFilePath.c_str(), "w");
	if (fpt == NULL) {
		printf("Error:unable to open Point file!\n");
		return false;
	}

	FILE* fdes = fopen(m_DesFilePath.c_str(), "w");
	if (fdes == NULL) {
		printf("Error:unable to open DesFilePath file!\n");
		return false;
	}
	FILE* fdes2;
	if ((fdes2 = fopen(m_DesFileTxt.c_str(), "w")) == NULL)
	{
		printf("can't open fdes2 file\n");
		return false;
	}

	map<int, PointXYZ>::iterator iter;
	int Point_id = 0;
	for (iter = m_AllPoint.begin(); iter != m_AllPoint.end(); iter++)
	{
		PointXYZ p = iter->second;
		if (p.n_Frame > 2 && p.n_Frame < 50 && p.VarXYZ[0] < 10 && p.VarXYZ[4] < 10 && p.VarXYZ[8] < 10)
		{
			//输出地图点文件 
			//fprintf(fpt, "%10d %10d %10d %13.4lf %13.4lf %13.4lf %13.4lf %13.4lf %13.4lf \n", iter->first, p.n_Frame, p.Recent_FI, p.XYZ[0], p.XYZ[1], p.XYZ[2], p.VarXYZ[0], p.VarXYZ[1], p.VarXYZ[2]);
			//提取相应的des
			vector<vector<uchar>> Des;
			for (int j = 0; j < p.n_Frame; j++)
			{
				int f_id = p.uvlist[j].Frame_id;
				int d_id = p.uvlist[j].Des_id;
				vector<uchar> des;
				des.resize(128);
				stringstream ss;
				ss << setfill('0') << setw(6) << f_id;
				string path = m_PicFolderPath + "\\SIFT\\" + ss.str() + ".sift";
				FILE* fd = fopen(path.c_str(), "rb");

				if (!fd)
				{
					printf("Fail to Open %s \n", path);
				}
				else
				{
					int Key_Num = 0;
					int Des_idx = 0;
					int Point_idx = 0;
					float u, v;
					fread(&Key_Num, sizeof(int), 1, fd);
					while (!feof(fd))
					{
						fread(&Des_idx, sizeof(int), 1, fd);
						fread(&Point_idx, sizeof(int), 1, fd);
						fread(&u, sizeof(float), 1, fd);
						fread(&v, sizeof(float), 1, fd);
						if (Des_idx == d_id && Point_idx == Point_id)
						{
							fread(&des[0], sizeof(uchar), 128, fd);
							break;
						}
						else
						{
							fseek(fd, 128L, SEEK_CUR);
						}
					}
					Des.push_back(des);
					des.clear();
				} 
				fclose(fd);
			}
			//SelectDesTmp(Des, fdes2);
			SelectDes(Des, fdes, fdes2);
		}
		Point_id++;
	}

	fclose(fpt);
	fclose(fdes);
	fclose(fdes2);
	return true;
}

bool CPointCloud::SelectPointTmp(FILE* fp, FILE* fdes, FILE* fp2, FILE* fdes2, int Point_id, PInfo *pinfo)
{
	PointXYZ P = m_AllPoint[Point_id];

	//if (pinfo->MMS_Ang[1] < 39.76 && (pinfo->Res < 1.85) && pinfo->MMS_Lat[1]<17.87 && pinfo->MMS_Dis[1] < 33.22)
    if (pinfo->MMS_Ang[1] < 39.76 && (pinfo->Res < 5) && pinfo->MMS_Lat[1] < 18 && pinfo->MMS_Dis[1] < 36)
	//if (pinfo->MMS_Ang[1] < 60 && (pinfo->Res < 9) && pinfo->MMS_Lat[1] < 36 && pinfo->MMS_Dis[1] < 72)
	//if (pinfo->MMS_Ang[1] < 90)
	{

	}
	else
	{
		return false;
	}
	
	double VarP[6] = { 0.0 };
	VarP[0] = P.VarXYZ[0];
	VarP[1] = P.VarXYZ[4];
	VarP[2] = P.VarXYZ[8];
	VarP[3] = P.VarXYZ[1];
	VarP[4] = P.VarXYZ[2];
	VarP[5] = P.VarXYZ[5];
	fwrite(P.XYZ, sizeof(double), 3, fp);
	fwrite(VarP, sizeof(double), 6, fp);
	fwrite(&P.n_Frame, sizeof(int), 1, fp);

	fprintf(fp2, "%13.4lf %13.4lf %13.4lf ", P.XYZ[0], P.XYZ[1], P.XYZ[2]);
	fprintf(fp2, "%13.4lf %13.4lf %13.4lf %13.4lf %13.4lf %13.4lf ", VarP[0], VarP[1], VarP[2], VarP[3], VarP[4], VarP[5]);
	fprintf(fp2, "%4d ", P.n_Frame);

	vector<vector<uchar>> vec_Des;

	for (int i = 0; i < P.n_Frame; i++)
	{
		fwrite(&m_ValidDesIndex, sizeof(int), 1, fp);
		fprintf(fp2, "%4d %4d", m_ValidPointIndex, P.uvlist[i].Frame_id);
		//寻找对应帧的描述子文件 打开并寻找对应的描述子
		int F_id = P.uvlist[i].Frame_id;
		int D_id = P.uvlist[i].Des_id;
		stringstream ss;
		ss << setfill('0') << setw(6) << F_id;
		string path = m_PicFolderPath + "\\SIFT\\" + ss.str() + ".sift";
		vector<uchar> Des;
		Des.resize(128);

		FILE* fd = fopen(path.c_str(), "rb");
		if (!fd)
		{
			printf("Fail to Open %s \n", path);
		}
		else
		{
			int Key_Num = 0;
			int Des_idx = 0;
			int Point_idx = 0;
			float u, v;
			fread(&Key_Num, sizeof(int), 1, fd);
			while (!feof(fd))
			{
				fread(&Des_idx, sizeof(int), 1, fd);
				fread(&Point_idx, sizeof(int), 1, fd);
				fread(&u, sizeof(float), 1, fd);
				fread(&v, sizeof(float), 1, fd);
				if (Des_idx == D_id && Point_idx == Point_id)
				{
					fread(&Des[0], sizeof(uchar), 128, fd);
					break;
				}
				else
				{
					fseek(fd, 128L, SEEK_CUR);
				}
			}
			vec_Des.push_back(Des);
		}

		fclose(fd);

		if (Des.size() != 128)
			printf("Search Descriptor wrong!\n");
		else
		{
			fwrite(&Des[0], sizeof(uchar), 128, fdes);
			for (int k = 0; k < 128; k++)
			{
				fprintf(fdes2, "%4u", Des[k]);
			}
			fprintf(fdes2, "\n");	
		}	
		m_ValidDesIndex++;
		
	}
	SelectDes(vec_Des, fdes, fdes2);
	vec_Des.clear();

	fprintf(fp2, "\n");
	m_ValidPointIndex++;

	fflush(fp);
	fflush(fp2);
	fflush(fdes);
	fflush(fdes2);
	return true;
}

void CPointCloud::SelectDes(vector<vector<uchar>> Des, FILE * fdes, FILE * fdes2)
{
	CKMeans kmeans;
	//????
	FILE* fp = fopen("D:\\KITTI\\04\\kmeans_h.txt", "a+");
	vector<double> H;
	double h;
	double h_sta[3];//聚类后各点距离的统计量
	double h_minmax[2];

	kmeans.SetInputDesDatabase(Des);
	if (Des.size() <= 5)
	{
		kmeans.setK(1);
	}
	else if (Des.size() > 5 && Des.size() <= 20)
	{
		kmeans.setK(2);
	}
	else
	{
		kmeans.setK(4);
	}
	kmeans.InitKCenter();
	kmeans.Cluster();	
	
	//kmeans.SaveFile(fdes, fdes2);

	fprintf(fp, "%d ", Des.size());
	fprintf(fp, "%d ", kmeans.m_k);

	for (int i = 0; i < kmeans.m_k; i++)
	{
		for (int j = 0; j < kmeans.m_grp_desDB[i].size(); j++)
		{
			h = kmeans.RealDistBetweenDes(kmeans.m_grp_desDB[i][j].des, kmeans.m_center[i]);
			H.push_back(h);
		}
	}
	GetMedMeanSTD(H, h_sta);
	GetMinMax(H, h_minmax);
	fprintf(fp, "%lf %lf %lf ", h_sta[0], h_sta[1], h_sta[2]);
	fprintf(fp, "%lf %lf ,", h_minmax[0], h_minmax[1]);

	//2019.03.08
	//同一特征点聚类中心的相互距离
	for (int i = 0, j = 1; j < kmeans.m_k; i++, j++)
	{
		h = kmeans.RealDistBetweenDes(kmeans.m_center[i], kmeans.m_center[j]);
		fprintf(fp, "%lf ", h);
	}
	fprintf(fp, ", ");
	//m_center_seq.push_back(kmeans.m_center);
	m_center_seq[m_ValidPointIndex] = kmeans.m_center;

	if (m_center_seq.size() > 10)
	{
		for (int j = 0; j < m_center_seq[m_ValidPointIndex].size(); j++)
		{
			for (int i = m_ValidPointIndex - 1; i >= m_ValidPointIndex - 10; i--)
			{
				for (int m = 0; m <m_center_seq[i].size(); m++)
				{
					h = kmeans.RealDistBetweenDes(m_center_seq[m_ValidPointIndex][j], m_center_seq[i][m]);
					fprintf(fp, "%lf ", h);
				}
			}
		}
	}
	else
	{
		for (int j = 0; j < m_center_seq[m_ValidPointIndex].size(); j++)
		{
			for (int i = m_ValidPointIndex - 1; i >= 0; i--)
			{
				for (int m = 0; m < m_center_seq[i].size(); m++)
				{
					h = kmeans.RealDistBetweenDes(m_center_seq[m_ValidPointIndex][j], m_center_seq[i][m]);
					fprintf(fp, "%lf ", h);
				}
			}
		}
	}
	fprintf(fp, ", ");
	//2019.03.08

	//fprintf(fp, "%lf %lf %lf ", h_sta[0], h_sta[1], h_sta[2]);
	//fprintf(fp, "%lf %lf ", h_minmax[0], h_minmax[1]);
	//for (int i = 0; i < H.size(); i++)
	//{
	//	fprintf(fp, "%lf ", H[i]);
	//}
	fprintf(fp, "\n");
	
	fflush(fdes2);
	fclose(fp);
}


void CPointCloud::SelectDesTmp(vector<vector<uchar>> Des, FILE * fdes)
{
	//int n = Des.size();
	//vector<uchar> Des_mean;
	//Des_mean.resize(128);
	//for (int i = 0; i < n; i++)
	//{
	//	for (int j = 0; j < 128; j++)
	//	{
	//		Des_mean[j] = Des_mean[j] + Des[i][j] / n;
	//	}		
	//}
	//for (int i = 0; i < n; i++)
	//{
	//	int dis = CalDesDis(Des_mean, Des[i]);
	//	fprintf(fdes, "%d\n", dis);
	//}


}

bool CPointCloud::SaveFile()
{
	FILE* fpt = fopen(m_PointFilePath.c_str(), "w");
	if (fpt == NULL) {
		printf("Error:unable to open Point file!\n");
		return false;
	}
	map<int, PointXYZ>::iterator iter;
	for (iter = m_AllPoint.begin(); iter != m_AllPoint.end(); iter++)
	{
		PointXYZ p = iter->second;
		if(p.VarXYZ[0]<10 && p.VarXYZ[1] < 10 && p.VarXYZ[2] < 10)
			fprintf(fpt, "%10d %10d %10d %13.4lf %13.4lf %13.4lf %13.4lf %13.4lf %13.4lf \n", iter->first, p.n_Frame, p.Recent_FI,p.XYZ[0], p.XYZ[1], p.XYZ[2], p.VarXYZ[0], p.VarXYZ[1], p.VarXYZ[2]);
	}
	fclose(fpt);
	
	return true;
}

//uv数组中为UL,V,UL-UR；
void CPointCloud::BacProject(double * uv, double * xyz)
{
	xyz[2] = m_bf / uv[2];
	xyz[0] = (uv[0] - m_cx)*xyz[2] / m_fx;
	xyz[1] = (uv[1] - m_cy)*xyz[2] / m_fy;
}

void CPointCloud::StereoProject(double * xyz, double * uv)
{
	uv[0] = m_fx * xyz[0] / xyz[2] + m_cx;
	uv[1] = m_fy * xyz[1] / xyz[2] + m_cy;
	uv[2] = m_bf / xyz[2];
}

void CPointCloud::ReadSiftInfo(int Frame_id, SIFTInfo * siftinfo)
{
	stringstream ss;
	ss << setfill('0') << setw(6) << Frame_id;
	string path = m_PicFolderPath + "\\SIFT\\" + ss.str() + ".sift";
	FILE* fd = fopen(path.c_str(), "rb");
	if (fd)
	{
		fread(&siftinfo->KeyNum, sizeof(int), 1, fd);
		siftinfo->Point_index.resize(siftinfo->KeyNum);
		siftinfo->Keys.resize(siftinfo->KeyNum);
		siftinfo->des.resize(siftinfo->KeyNum * 128);

		for (int i = 0; i < siftinfo->KeyNum; i++)
		{
			int Des_idx;
			fread(&Des_idx, sizeof(int), 1, fd);
			fread(&siftinfo->Point_index[i], sizeof(int), 1, fd);
			fread(&siftinfo->Keys[i].x, sizeof(float), 1, fd);
			fread(&siftinfo->Keys[i].y, sizeof(float), 1, fd);
			fread(&siftinfo->des[i*128], sizeof(uchar), 128, fd);
		}
	}
	fclose(fd);
}

void CPointCloud::GetMedMeanSTD(vector<double> V, double* MMS)
{
	double sum = 0.0;
	double var = 0.0;
	int n = V.size();
	sort(V.begin(), V.end());
	if (n % 2)  
		MMS[0] = V[floor(n / 2)]; 
	else 
		MMS[0] = (V[n / 2] + V[n / 2 - 1]) / 2;

	for (int i = 0; i < n; i++) sum = sum + V[i];

	MMS[1] = sum / n;
	for (int i = 0; i < n; i++)
	{
		var = var + (V[i] - MMS[1])*(V[i] - MMS[1]);
	}
	MMS[2] = sqrt(var / n);
}

void CPointCloud::GetMinMax(vector<double> v, double * MM)
{
	vector<double>::iterator min,max;

	min = min_element(begin(v), end(v));
	max = max_element(begin(v), end(v));
	MM[0] = *min;
	MM[1] = *max;
}

int CPointCloud::CalDesDis(vector<uchar> des1, vector<uchar> des2)
{
	int Sum = 0;
	for (int i = 0; i < 128; i++)
	{
		Sum = Sum + (des1[i] - des2[i])*(des1[i] - des2[i]);
	}
	return sqrt(Sum);
}

void CPointCloud::CalBackProjectError()
{
	for (int i = 0; i < CurPoint.uvlist.size(); i++)
	{
		if (CurPoint.uvlist[i].ValidFlag == 0)
			continue;

		double P_tmp[3], P_tmpX[9], p[3], uvd[3];

		int F_id = CurPoint.uvlist[i].Frame_id;
		MatrixSubtraction(3, 1, CurPoint.XYZ, m_vGroundTruth[F_id].Pos, P_tmp);
		MatrixSkewSymmetric(P_tmp, P_tmpX);
		MatrixMultiply(3, 3, m_vGroundTruth[F_id].Rwc, 3, 1, P_tmp, p);
		StereoProject(p, uvd);
		MatrixSubtraction(3, 1, CurPoint.uvlist[i].UV, uvd, CurPoint.uvlist[i].BackProjectUVError);
		CurPoint.uvlist[i].MeanBP_UVError= CurPoint.uvlist[i].MeanBP_UVError = (fabs(CurPoint.uvlist[i].BackProjectUVError[0]) +
			fabs(CurPoint.uvlist[i].BackProjectUVError[1]) +
			fabs(CurPoint.uvlist[i].BackProjectUVError[2])) / 3;
		if (CurPoint.uvlist[i].MeanBP_UVError > m_BackProjectUVError)
			m_BAState = 0;
	}
}


