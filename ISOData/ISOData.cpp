// ISOData.cpp : �������̨Ӧ�ó������ڵ㡣
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

#define iniClusters 5  //��ʼ��۵ĸ���
typedef int int32;
using namespace std;

//����6��ʹ�õĲ���
struct Args
{
	int expClusters;   //�����õ��ľ�����
	int thetaN;        //����������������
	int maxIts;        //����������
	int combL;         //ÿ�ε�������ϲ������������
	double thetaS;     //��׼ƫ�����
	double thetaC;     //�ϲ�����
}args;

//�����ά�㣬��������Ƕ�ά����������Ȼ�����ƹ㵽��ά
struct Point
{
	double x, y;
};

//��Ҫ�ϲ���������۵���Ϣ������������۵�id�;���
struct MergeInfo
{
	int u, v;
	double d;    //���u���������v���ĵľ���
};

//����ȽϺ���
bool cmp(MergeInfo a, MergeInfo b)
{
	return a.d < b.d;
}

//��������֮�����
double dist(Point A, Point B)
{
	return sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y));
}

struct Cluster
{
	int nSamples;          //������ĸ���
	double avgDist;        //�����㵽�������ĵ�ƽ������
	Point center;          //��������
	Point sigma;           //���������ĵı�׼��
	vector<Point *> data;  //���������

	//����þ�������ģ�������ľ�ֵ
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

	//������������㵽�þ������ĵ�ƽ������
	void calDist()
	{
		avgDist = 0;
		for(int i = 0; i < nSamples; i++)
			avgDist += dist(*(data.at(i)), center);
		avgDist /= nSamples;
	}

	//�������������ĵı�׼��
	void calStErr()
	{
		assert(nSamples == data.size());
		double attr1 = 0;
		double attr2 = 0;        //����������ά��
		for(int i = 0; i < nSamples; i++)
		{
			attr1 += (data.at(i)->x - center.x) * (data.at(i)->x - center.x);
			attr2 += (data.at(i)->y - center.y) * (data.at(i)->y - center.y);
		}
		sigma.x = sqrt(attr1 / nSamples);
		sigma.y = sqrt(attr2 / nSamples);
	}
};

//��ȡ����
void getData(Point p[], int n)
{
	cout << "getting data..." << endl;
	for(int i = 0; i < n; i++)
		scanf("%lf %lf", &p[i].x, &p[i].y);
	cout << "get data done!" << endl;
}

//���ò�����ֵ
void setArgs()
{
	args.expClusters = 5;
	args.thetaN = 3;
	args.maxIts = 10000;
	args.combL = 10;
	args.thetaS = 3;
	args.thetaC = 0.001;
}

//Ѱ�ҵ�t���������������Ķ�Ӧ��id
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

//���ַ�Ѱ�Ҿ���պ�С��thetaC��������۵�index
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
		cout << "��" << i + 1 << "��������:" << endl;
		for(int j = 0; j < c.at(i).data.size(); j++)
			cout << "(" << c[i].data[j]->x << "," << c[i].data[j]->y << ")  ";
		cout << endl;
		cout << endl;
	}
}

void ISOData(Point p[], int n)
{
	cout << "ISOData is processing......." << endl;
	vector<Cluster> c;              //ÿ����۵�����
	const double split = 0.5;       //���ѳ���(0,1]
	int nClusters = iniClusters;    //��ʼ����۸���

	//��ʼ��nClusters���࣬�����������
	for(int i = 0; i < nClusters; i++)
	{
		Cluster t;
		t.center = p[i];
		t.nSamples = 0;
		t.avgDist = 0;
		c.push_back(t);
	}

	int iter = 0;
	bool isLess = false;            //��־�Ƿ��������Ŀ����thetaN
	while(1)
	{
		//�����ÿһ������
		for(int i = 0; i < nClusters; i++)
		{
			c.at(i).nSamples = 0;
			c.at(i).data.clear();
		}

		//�������������ֵ���������������������
		for(int i = 0; i < n; i++)
		{
			int idx = FindIdx(c, p[i]);
			c.at(idx).data.push_back(&p[i]);
			c.at(idx).nSamples++;
		}

		int k = 0;                   //��¼������Ŀ����thetaN�����index
		for(int i = 0; i < nClusters; i++)
		{
			if(c.at(i).data.size() < args.thetaN)
			{
				isLess = true;       //˵�����������٣�����Ӧ��ɾ��
				k = i;
				break;
			}
		}

		//��������������ĿС��thetaN
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

		//���¼����ֵ��������������ĵ�ƽ������
		for(int i = 0; i < nClusters; i++)
		{
			c.at(i).calMean();
			c.at(i).calDist();
		}

		//�����ܵ�ƽ������
		double totalAvgDist = 0;
		for(int i = 0; i < nClusters; i++)
			totalAvgDist += c.at(i).avgDist * c.at(i).nSamples;
		totalAvgDist /= n;

		if(iter >= args.maxIts) break;

		//���Ѳ���
		if(nClusters <= args.expClusters / 2)
		{
			vector<double> maxsigma;
			for(int i = 0; i < nClusters; i++)
			{
				//�������ı�׼ƫ��
				c.at(i).calStErr();
				//��������׼���������
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
						Cluster newCtr;     //�µľ�������
						//��ȡ�µ�����
						newCtr.center.x = c.at(i).center.x - split * c.at(i).sigma.x;
						newCtr.center.y = c.at(i).center.y - split * c.at(i).sigma.y;
						c.push_back(newCtr);
						//�ı��ϵ�����
						c.at(i).center.x = c.at(i).center.x + split * c.at(i).sigma.x;
						c.at(i).center.y = c.at(i).center.y + split * c.at(i).sigma.y;
						break;
					}
				}
			}
		}

		//�ϲ�����
		if(nClusters >= 2 * args.expClusters || (iter & 1) == 0)
		{
			int size = nClusters * (nClusters - 1);
			//��Ҫ�ϲ��ľ������
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
			//��������
			sort(info, info + cnt, cmp);
			//�ҳ�info�����о���պ�С��thetaC��index����ôindex��С�ĸ�Ӧ�úϲ�
			int iPos = FindPos(info, cnt, args.thetaC);

			//����ָʾ��λ�õ��������Ƿ��Ѿ��ϲ�
			bool *flag = new bool[nClusters];
			memset(flag, false, sizeof(bool) * nClusters);
			//���ڱ�Ǹ�λ�õ��������Ƿ��Ѿ��ϲ�ɾ��
			bool *del = new bool[nClusters];
			memset(del, false, sizeof(bool) * nClusters);
			//��¼�ϲ��Ĵ���
			int nTimes = 0;

			for(int i = 0; i <= iPos; i++)
			{
				int u = info[i].u;
				int v = info[i].v;
				//ȷ��ͬһ�����ֻ�ϲ�һ��
				if(!flag[u] && !flag[v])
				{
					nTimes++;
					//���һ�ε����кϲ���������combL����ֹͣ�ϲ�
					if(nTimes > args.combL) break;
					//����Ŀ�ٵ������ϲ�����Ŀ���������
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

			//ɾ���ϲ���ľ���
			vector<Cluster>::iterator id = c.begin();
			for(int i = 0; i < nClusters; i++)
			{
				if(del[i])
					id = c.erase(id);
				else
					id++;
			}

			//�ϲ����ٴξ�ɾ�����ٸ�
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

