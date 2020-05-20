#pragma once
#ifndef PCOPTION_H
#define PCOPTION_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <BaseSDC.h>
#include <BaseMath.h>
#include <BaseCmnFunc.h>
#include <BaseMatrix.h>

using namespace std;
using namespace basetk;
class CPCOPTION
{
public:
	CPCOPTION();
	virtual ~CPCOPTION();
	bool readOptFile(string file);
	

	float m_fx;
	float m_fy;
	float m_cx;
	float m_cy;
	float m_bf;
	int   m_StartFid;
	int   m_SuffixFlag;
	double m_VarPix;
	double m_BackProjectUVError;
	string m_DesFilePath;
	string m_PointFilePath;
	string m_UVFilePath;
	string m_PointInfoPath;
	string m_TimesPath;
	string m_GTFilePath;
	string m_PicLFolderPath;
	string m_PicRFolderPath;
	string m_PICFolder;
private:
	int  DecodeEachOpt(const char* oneline); // 0:正常,1:数据读完,2:没有相应Opt
	

};



#endif

