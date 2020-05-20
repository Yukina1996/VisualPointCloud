#include "SIFT.h"
#include "PCOption.h"


#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include <iostream>
#include <time.h>

#define CLOCKS_PER_SEC ((clock_t)1000)

using namespace basetk;


//步骤一：从序列影像和真值位姿中增量构建特征点的三维地图数据

//步骤二：将构建出的每个地图点采用不同的像素观测值进行平差

int main()
{
	clock_t start = clock();
	string Optfile = "D:\\PointCloud\\PointCloud\\Option.txt";

	CPointCloud PC;	
	PC.runPointCloud(Optfile);

	printf("over!\n");
	clock_t end = (clock() - start);
	cout << "time comsumption is " << end << endl;
	getchar();
	return 1;
}




