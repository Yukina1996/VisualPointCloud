#pragma once
#ifndef SIFT_H
#define SIFT_H

#include <SiftGPU\SiftGPU.h>
#include <map>
#include "PCOption.h"
#include "base.h"
#include <direct.h>
#include <CPose3D.h>



class CPointCloud
{
public:
	CPointCloud();
	virtual ~CPointCloud();
	//步骤一：从序列影像和真值位姿中增量构建特征点的三维地图数据
	bool runPointCloud(string optfile, CPCOPTION* pPCOpt = NULL);
private:
	//step 1:
	bool Init(CPCOPTION* pPCOpt);
	bool readGTFile(string path);
	//step 2:
	bool RunOneLRSIFT();
	bool InitPoint();
	void RunSeqMatch();
	double EssentialVerify(int Frameidx, double Cur_uv[], double Match_uv[]);
	void UpdateInfo();
	//step 3:
	int BAPoint(int Point_id);
	int RunOneBA();
	bool DealAllPoint();
    //step 6:
	void PointInfo(FILE* fout_uv, FILE* fout_pi, int Point_id, PInfo *pinfo);
	//step 7:
	bool SelectPoint();
	bool SelectPointTmp(FILE* fp, FILE* fdes, FILE* fp2, FILE* fdes2, int Point_id, PInfo *pinfo);

	void SelectDes(vector<vector<uchar>> Des,FILE *fdes, FILE *fdes2);
	void SelectDesTmp(vector<vector<uchar>> Des, FILE *fdes);

	bool SaveFile();

	void BacProject(double* uv, double* xyz);
	void StereoProject(double* xyz, double* uv);
	void ReadSiftInfo(int Frame_id, SIFTInfo *siftinfo);
	void GetMedMeanSTD(vector<double> V, double* MMS);
	void GetMinMax(vector<double> v, double* MM);
	int CalDesDis(vector<uchar>des1, vector<uchar>des2);
	void CalBackProjectError();

public:
	vector<string> m_vstrImageLeft;
	vector<string> m_vstrImageRight;
	vector<GPSTIME> m_vTimestamps;
	vector<GroundTruth> m_vGroundTruth;
	PointXYZ CurPoint;

	//vector<vector<vector<uchar>>> m_center_seq;
	map<int, vector<vector<uchar>>> m_center_seq;

private:
	CPCOPTION   m_PCOpt;
	float m_fx;
	float m_fy;
	float m_cx;
	float m_cy;
	float m_bf;
	double m_BackProjectUVError;
	int m_BAState;

	
	double m_VarPix;
	string m_PicFolderPath;
	string m_DesFilePath;
	string m_PointFilePath;
	string m_UVFilePath;
	string m_PointFileTxt;
	string m_PointInfoPath;
	string m_DesFileTxt;

	int m_frameIndex;
	int m_PointIndex;
	int m_DesIndex;

	int m_ValidDesIndex;
	int m_ValidPointIndex;

	double m_CurFrameAtt[9];  //Rwc
	double m_CurFramePos[3];
	SIFTInfo m_Cursiftinfo;
	map<int,PointXYZ> m_AllPoint;//地图点编号从0开始
};





#endif // !SIFT_H
