// ISOData.cpp : 定义控制台应用程序的入口点。
//
#include <iostream>
using namespace std;
#include "stdafx.h"

#include <iostream>
#include <string.h>
#include <algorithm>
#include <stdio.h>
#include <vector>
#include <assert.h>
#include <math.h>

#define iniClusters 5  //初始类聚的个数
typedef int int32;
using namespace std;

//定义6个使用的参数
struct Args
{
	int expClusters;   //期望得到的聚类数
	int thetaN;        //聚类中最少样本数
	int maxIts;        //最大迭代次数
	int combL;         //每次迭代允许合并的最大聚类对数
	double thetaS;     //标准偏差参数
	double thetaC;     //合并参数
}args;

//定义二维点，这里假设是二维的特征，当然可以推广到多维
struct Point
{
	double x, y;
};

//需要合并的两个类聚的信息，包括两个类聚的id和距离
struct MergeInfo
{
	int u, v;
	double d;    //类聚u中心与类聚v中心的距离
};

//定义比较函数
bool cmp(MergeInfo a, MergeInfo b)
{
	return a.d < b.d;
}

//计算两点之间距离
double dist(Point A, Point B)
{
	return sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y));
}

struct Cluster
{
	int nSamples;          //样本点的个数
	double avgDist;        //样本点到样本中心的平均距离
	Point center;          //样本中心
	Point sigma;           //样本与中心的标准差
	vector<Point *> data;  //聚类的数据

	//计算该聚类的中心，即该类的均值
	void calMean()
	{
		assert(nSamples == data.size());
		for(int i = 0; i < nSamples; i++)
		{
			center.x += data.at(i)->x;
			center.y += data.at(i)->y;
		}
		center.x /= nSamples;
		center.y /= nSamples;
	}

	//计算该类样本点到该聚类中心得平均距离
	void calDist()
	{
		avgDist = 0;
		for(int i = 0; i < nSamples; i++)
			avgDist += dist(*(data.at(i)), center);
		avgDist /= nSamples;
	}

	//计算样本与中心的标准差
	void calStErr()
	{
		assert(nSamples == data.size());
		double attr1 = 0;
		double attr2 = 0;        //样本的两个维度
		for(int i = 0; i < nSamples; i++)
		{
			attr1 += (data.at(i)->x - center.x) * (data.at(i)->x - center.x);
			attr2 += (data.at(i)->y - center.y) * (data.at(i)->y - center.y);
		}
		sigma.x = sqrt(attr1 / nSamples);
		sigma.y = sqrt(attr2 / nSamples);
	}
};

//获取数据
void getData(Point p[], int n)
{
	cout << "getting data..." << endl;
	for(int i = 0; i < n; i++)
		scanf("%lf %lf", &p[i].x, &p[i].y);
	cout << "get data done!" << endl;
}

//设置参数的值
void setArgs()
{
	args.expClusters = 5;
	args.thetaN = 3;
	args.maxIts = 10000;
	args.combL = 10;
	args.thetaS = 3;
	args.thetaC = 0.001;
}

//寻找点t距离最近的类的中心对应的id
int FindIdx(vector<Cluster> &c, Point &t)
{
	int nClusters = c.size();
	assert(nClusters >= 1);
	double ans = dist(c.at(0).center, t);
	int idx = 0;
	for(int i = 1; i < nClusters; i++)
	{
		double tmp = dist(c.at(i).center, t);
		if(ans > tmp)
		{
			idx = i;
			ans = tmp;
		}
	}
	return idx;
}

//二分法寻找距离刚好小于thetaC的两个类聚的index
int FindPos(MergeInfo *info, int n, double thetaC)
{
	int l = 0;
	int r = n - 1;
	while(l <= r)
	{
		int mid = (l + r) >> 1;
		if(info[mid].d < thetaC)
		{
			l = mid + 1;
			if(l < n && info[l].d >= thetaC)
				return mid;
		}
		else
		{
			r = mid - 1;
			if(r >= 0 && info[r].d < thetaC)
				return r;
		}
	}
	if(info[n - 1].d < thetaC)
		return n - 1;
	else
		return -1;
}

void Print(const vector<Cluster> c)
{
	int n = c.size();
	for(int i = 0; i < n; i++)
	{
		cout << "------------------------------------" << endl;
		cout << "第" << i + 1 << "个聚类是:" << endl;
		for(int j = 0; j < c.at(i).data.size(); j++)
			cout << "(" << c[i].data[j]->x << "," << c[i].data[j]->y << ")  ";
		cout << endl;
		cout << endl;
	}
}

void ISOData(Point p[], int n)
{
	cout << "ISOData is processing......." << endl;
	vector<Cluster> c;              //每个类聚的数据
	const double split = 0.5;       //分裂常数(0,1]
	int nClusters = iniClusters;    //初始化类聚个数

	//初始化nClusters个类，设置相关数据
	for(int i = 0; i < nClusters; i++)
	{
		Cluster t;
		t.center = p[i];
		t.nSamples = 0;
		t.avgDist = 0;
		c.push_back(t);
	}

	int iter = 0;
	bool isLess = false;            //标志是否有类的数目低于thetaN
	while(1)
	{
		//先清空每一个聚类
		for(int i = 0; i < nClusters; i++)
		{
			c.at(i).nSamples = 0;
			c.at(i).data.clear();
		}

		//将所有样本划分到距离类聚中心最近的类中
		for(int i = 0; i < n; i++)
		{
			int idx = FindIdx(c, p[i]);
			c.at(idx).data.push_back(&p[i]);
			c.at(idx).nSamples++;
		}

		int k = 0;                   //记录样本数目低于thetaN的类的index
		for(int i = 0; i < nClusters; i++)
		{
			if(c.at(i).data.size() < args.thetaN)
			{
				isLess = true;       //说明样本数过少，该类应该删除
				k = i;
				break;
			}
		}

		//如果有类的样本数目小于thetaN
		if(isLess)
		{
			nClusters--;
			Cluster t = c.at(k);
			vector<Cluster>::iterator pos = c.begin() + k;
			c.erase(pos);
			assert(nClusters == c.size());
			for(int i = 0; i < t.data.size(); i++)
			{
				int idx = FindIdx(c, *(t.data.at(i)));
				c.at(idx).data.push_back(t.data.at(i));
				c.at(idx).nSamples++;
			}
			isLess = false;
		}

		//重新计算均值和样本到类聚中心的平均距离
		for(int i = 0; i < nClusters; i++)
		{
			c.at(i).calMean();
			c.at(i).calDist();
		}

		//计算总的平均距离
		double totalAvgDist = 0;
		for(int i = 0; i < nClusters; i++)
			totalAvgDist += c.at(i).avgDist * c.at(i).nSamples;
		totalAvgDist /= n;

		if(iter >= args.maxIts) break;

		//分裂操作
		if(nClusters <= args.expClusters / 2)
		{
			vector<double> maxsigma;
			for(int i = 0; i < nClusters; i++)
			{
				//计算该类的标准偏差
				c.at(i).calStErr();
				//计算该类标准差的最大分量
				double mt = c.at(i).sigma.x > c.at(i).sigma.y? c.at(i).sigma.x : c.at(i).sigma.y;
				maxsigma.push_back(mt);
			}
			for(int i = 0; i < nClusters; i++)
			{
				if(maxsigma.at(i) > args.thetaS)
				{
					if((c.at(i).avgDist > totalAvgDist && c.at(i).nSamples > 2 * (args.thetaN + 1)) || (nClusters < args.expClusters / 2))
					{
						nClusters++;
						Cluster newCtr;     //新的聚类中心
						//获取新的中心
						newCtr.center.x = c.at(i).center.x - split * c.at(i).sigma.x;
						newCtr.center.y = c.at(i).center.y - split * c.at(i).sigma.y;
						c.push_back(newCtr);
						//改变老的中心
						c.at(i).center.x = c.at(i).center.x + split * c.at(i).sigma.x;
						c.at(i).center.y = c.at(i).center.y + split * c.at(i).sigma.y;
						break;
					}
				}
			}
		}

		//合并操作
		if(nClusters >= 2 * args.expClusters || (iter & 1) == 0)
		{
			int size = nClusters * (nClusters - 1);
			//需要合并的聚类个数
			int cnt = 0;
			MergeInfo *info = new MergeInfo[size];
			for(int i = 0; i < nClusters; i++)
			{
				for(int j = i + 1; j < nClusters; j++)
				{
					info[cnt].u = i;
					info[cnt].v = j;
					info[cnt].d = dist(c.at(i).center, c.at(j).center);
					cnt++;
				}
			}
			//进行排序
			sort(info, info + cnt, cmp);
			//找出info数组中距离刚好小于thetaC的index，那么index更小的更应该合并
			int iPos = FindPos(info, cnt, args.thetaC);

			//用于指示该位置的样本点是否已经合并
			bool *flag = new bool[nClusters];
			memset(flag, false, sizeof(bool) * nClusters);
			//用于标记该位置的样本点是否已经合并删除
			bool *del = new bool[nClusters];
			memset(del, false, sizeof(bool) * nClusters);
			//记录合并的次数
			int nTimes = 0;

			for(int i = 0; i <= iPos; i++)
			{
				int u = info[i].u;
				int v = info[i].v;
				//确保同一个类聚只合并一次
				if(!flag[u] && !flag[v])
				{
					nTimes++;
					//如果一次迭代中合并对数多于combL，则停止合并
					if(nTimes > args.combL) break;
					//将数目少的样本合并到数目多的样本中
					if(c.at(u).nSamples < c.at(v).nSamples)
					{
						del[u] = true;
						Cluster t = c.at(u);
						assert(t.nSamples == t.data.size());
						for(int j = 0; j < t.nSamples; j++)
							c.at(v).data.push_back(t.data.at(j));
						c.at(v).center.x = c.at(v).center.x * c.at(v).nSamples + t.nSamples * t.center.x;
						c.at(v).center.y = c.at(v).center.y * c.at(v).nSamples + t.nSamples * t.center.y;
						c.at(v).nSamples += t.nSamples;
						c.at(v).center.x /= c.at(v).nSamples;
						c.at(v).center.y /= c.at(v).nSamples;
					}
					else
					{
						del[v] = true;
						Cluster t = c.at(v);
						assert(t.nSamples == t.data.size());
						for(int j = 0; j < t.nSamples; j++)
							c.at(u).data.push_back(t.data.at(j));
						c.at(u).center.x = c.at(u).center.x * c.at(u).nSamples + t.nSamples * t.center.x;
						c.at(u).center.y = c.at(u).center.y * c.at(u).nSamples + t.nSamples * t.center.y;
						c.at(u).nSamples += t.nSamples;
						c.at(u).center.x /= c.at(u).nSamples;
						c.at(u).center.y /= c.at(u).nSamples;
					}
				}
			}

			//删除合并后的聚类
			vector<Cluster>::iterator id = c.begin();
			for(int i = 0; i < nClusters; i++)
			{
				if(del[i])
					id = c.erase(id);
				else
					id++;
			}

			//合并多少次就删除多少个
			nClusters -= nTimes;
			assert(nClusters == c.size());
			delete[] info;
			delete[] flag;
			delete[] del;
			info = NULL;
			flag = NULL;
			del = NULL;
		}

		if(iter >= args.maxIts) break;
		iter++;
	}
	assert(nClusters == c.size());
	Print(c);
}

int _tmain(int argc, _TCHAR* argv[])
{
	int n;
	scanf("%d", &n);
	Point *p = new Point[n];

	getData(p, n);
	setArgs();
	ISOData(p, n);

	delete[] p;
	p = NULL;
	return 0;
}

