#ifndef PARAMETER_H
#define PARAMETER_H

// this function is used to input the parameter

#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;

//void negread(double & r0,double & rmax, int & imax, int & jmax,int & nulod ,double  & expan, vector<int> & npxy,
//				vector<vector<int> > & nelprt, int & ishape, int & ntype);
// ostream & operator<< ( vector<int> & ivec);

class parameter
{
    public:
        parameter();
        virtual ~parameter();

        void input();
        int check1();
        void echo();
        void negread();
        int nodexy();
        friend ostream & operator<< ( ostream& cout,vector<int> & ivec);
        void gaussq();
        int check2();

    protected:

    public:
        int npoin;
        int nelem;
        int nvfix;
        int ntype;
        int nnode;
        int nmats;
        int ngaus;
        int nalgo;
        int nincs;
        int nstre;
        int ncrit;
        int mcmky;
        int ndofn;
        int mfron;  //?  在什么地方定义的这个变量，在check2中检查了这个量


        int ishape, imax, jmax, nulod;
        double r0,rmax, expan, prstr, ptstr, hstre, vstre, zstre;
        //lack nnu ntype

        vector<int> npxy;	//？　与单元数量有关
        vector<vector<int> > nelprt ;	//存放x,y 的单元编号
        vector<int> matno ;	// ?    在dimen, input , zero 函数中均有调用
        vector<vector<int> > lnods ;        //?
        vector<vector<double> >  coord;   // nodal coordinates
        vector<int> nofix;      //?
        vector<vector<double> >  presc;  //?
        vector<int> iffix;      //?
        int nprop;
        vector<vector<double> >  props;     //? 存放物性参数
        vector<double>  posgp;         // gauss积分常数
        vector<double>  weigp;          //gauss 积分权重
        vector<int> ndfro ; //？ 在check2()函数中初始化；

};



#endif // PARAMETER_H
