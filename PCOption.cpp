#include "PCOption.h"

CPCOPTION::CPCOPTION()
{
	m_fx = 0.0;
	m_fy = 0.0;
	m_cx = 0.0;
	m_cy = 0.0;
	m_bf = 0.0;
	m_StartFid = 0;
	m_SuffixFlag = 0;
	m_VarPix = 0.0;
	m_DesFilePath = "";
	m_PointFilePath = "";
	m_UVFilePath = "";
	m_TimesPath = "";
	m_GTFilePath = "";

	m_PICFolder = "";
	m_PicLFolderPath = "";
	m_PicRFolderPath = "";
	
}

CPCOPTION::~CPCOPTION()
{
	
}

bool CPCOPTION::readOptFile(string file)
{
	FILE* fopt = fopen(file.c_str(), "rt");
	if (fopt == NULL) {
		printf("Error:unable to open opt file!\n"); 
		return false;
	}

	char oneline[1024] = { '\0' };
	while (fgets(oneline, 1024, fopt))
	{
		// ASCII 10 对应为换行符 LF 
		if (oneline[0] == 10 || oneline[0] == '#')
			continue;

		DecodeEachOpt(oneline);
	}

	fclose(fopt);

	return true;
}



int CPCOPTION::DecodeEachOpt(const char * oneline)
{
	if (strstr(oneline, "END"))  return 1;

	string text, head;
	int i = 0;
	if (CheckOptFileOneline(oneline, head, text) == false)
		return 2;

	char ch[1024] = { '\0' };

	if (head == string("Pic_Folder"))
	{
		sscanf(text.c_str(), "%s", ch);
		CheckFilePath(ch);
		m_PICFolder = string(ch);
		m_TimesPath = m_PICFolder + "\\times.txt";
		m_GTFilePath = m_PICFolder + "\\GroundTruth.txt";
		m_PicLFolderPath = m_PICFolder + "\\image_0\\";
		m_PicRFolderPath = m_PICFolder + "\\image_1\\";

		m_DesFilePath = m_PICFolder + "\\Descriptor.bin";
		m_PointFilePath = m_PICFolder + "\\Point.bin";
		m_UVFilePath = m_PICFolder + "\\UVlist.txt";
		m_PointInfoPath = m_PICFolder + "\\PointInfo.txt";

	}
	else if (head == string("Fx"))
	{
		sscanf(text.c_str(), "%f", &m_fx);
	}
	else if (head == string("Fy"))
	{
		sscanf(text.c_str(), "%f", &m_fy);
	}
	else if (head == string("Cx"))
	{
		sscanf(text.c_str(), "%f", &m_cx);
	}
	else if (head == string("Cy"))
	{
		sscanf(text.c_str(), "%f", &m_cy);
	}
	else if (head == string("Bf"))
	{
		sscanf(text.c_str(), "%f", &m_bf);
	}
	else if (head == string("VarPix"))
	{
		sscanf(text.c_str(), "%lf", &m_VarPix);
	}
	else if (head == string("BackProjectUVError"))
	{
		sscanf(text.c_str(), "%lf", &m_BackProjectUVError);
	}
	else if (head == string("Start_Fid"))
	{
		sscanf(text.c_str(), "%d", &m_StartFid);
	}
	else if (head == string("Suffix"))
	{
		sscanf(text.c_str(), "%d", &m_SuffixFlag);
	}
	return 0;
}


