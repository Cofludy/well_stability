#include "parameter.h"

parameter::parameter()
{
    //to initial the parameter

    npoin=3696 ;	//计算用节点总数
	nelem=1200;	    //计算用单元总数
	nvfix=96;	    //边界约束节点数
	ntype=2;	    //计算约束类型　＝１:平面应变；　＝２:平面应力；　
					//				＝３:轴对称；　＝４：完全平面应变
	nnode=8;	    //每个单元的节点数　＝4:线性等参四边形单元：
					//		　			＝８:serendipity 二次等参单元
	nmats=1;	    //计算用材料总数
	ngaus=2;	    // 高斯积分法制的选取，＝２:２点高斯积分法则；　＝３：３点高斯积分法则
	nalgo=1;	    //控制非线性问题的求解参数
 					//	＝１：初刚度法；　　＝２：线性刚度法
					//	＝３：联合算法，只对载荷增量的第一次迭代重新计算刚度矩阵
					//	＝４：联合算法，只对载荷增量的第二次迭代重新计算刚度矩阵（第一次除外）
	nincs=74;	    //达到最后增量的载荷总次数
	nstre=3;	    //问题的独立应力分量数，＝３平面问题；　＝４：　轴对称问题；　=5: 完全平面应变；
	ncrit=3;	    //屈服准则控制参数：
					//	=1:Trescr 屈服准则；　＝２：von Mises屈服准则
					//	=3: M-C屈服准则；　=4:Druker Prager 屈服准则；
	mcmky=1;	    // 对Ｍ－Ｃ模型，＝１: 考虑中间应力；　＝２：不考虑中间应力
					//对D-P模型：
					//＝１：　D-P模型和 M-C 模型是内接圆
					//＝２：　D-P模型和 M-C 模型是外接圆
					//＝３：　D-P模型和 M-C 模型是中接圆

    ishape = 1;
	imax=48 ;	    //周向单元数；
	jmax=25;	    //径向单元数
	nulod=5;	    //?
	r0=15.0;		//井孔半径
	expan=1.15;		//径向单元尺度的放大因子
	rmax=r0*11.0;	// ?　网格半径最大值
	prstr=0; ptstr=0;		// ?
	hstre=2000.0; vstre= 3000;  zstre=4000;	  //原地应力

	ndofn= 2;	    //？=2 for popular plane strain and
                        //  =3 for complete plane strain

	//vector<int> npxy(2,0);	//？　与单元数量有关
	for (int i=0;i!=2;++i)
    {
        npxy.push_back(0);
    }

    //vector<vector<int> > nelprt(2,vector<int> (50,0));	//存放x,y 的单元编号
    vector<int> tempa(50,0);
    for(int i=0;i!=2;++i)
    {
        nelprt.push_back(tempa);
    }

    //vector<int> matno(1300,0);	// ?    在dimen, input , zero 函数中均有调用
    for (int i=0;i!=1300;++i)
    {
        matno.push_back(0);
    }

    //vector<vector<int> > lnods(1300,vector<int> (9,0));
    vector<int> tempb(9,0);
    for(int i=0;i!=1300;++i)
    {
        lnods.push_back(tempb);
    }

    //vector<vector<double> >  coord(4900,2);
    vector<double> tempc(2,0);
    for(int i=0;i!=4900;++i)
    {
        coord.push_back(tempc);
    }

    //vector<int> nofix(288);
    for(int i=0;i!=288;++i)
    {
        nofix.push_back(0);
    }
    //vector<vector<double> >  presc(288,vector<double> (3));
    vector<double> tempd(3,0);
    for(int i=0;i!=288;i++)
    {
        presc.push_back(tempd);
    }
    //vector<int> iffix (14700)
    for(int i=0;i!=14700;++i)
    {
        iffix.push_back(0);
    }

    nprop = 7;        // 是props 的纵坐标
    //vector<vector<double>> props(5,vector<dopuble> (7))
    vector<double>  tempe(7,0);
    for(int i=0;i!=5;++i)
    {
        props.push_back(tempe);
    }

}

parameter::~parameter()
{
    //to delete the class
    cout<<" the parameter class was deleted! "<<endl;

}

void parameter::input( )
// to get the parameter
//输入给定的几何尺寸，边界条件，和材料性能数据
//为了方便传递一些参数放在向量中
{
    ifstream fcin;
	fcin.open("plas4825.in");
	if(!fcin)
	{
		cout<<"error:cannot open plas4852.in "<<endl;
		exit (1);
	}

	string title;
	getline(fcin,title);

    //read from the file "plas4825.in"
	fcin>>npoin>>nelem>>nvfix>>ntype>>nnode>>nmats>>ngaus>>nalgo>>nincs>>nstre>>ncrit>>mcmky;

    int ndofn(2);	//？=2 for popular plane strain and
					//  =3 for complete plane strain
	int nevab(0);	//? 单位矩阵的维数
	nevab = ndofn*nnode;

	int nstr1 = nstre + 1 ;	//独立应力分量数 + 1
	if (ntype == 3 )
		nstr1=nstre;

	int ntotv = npoin*ndofn;	//总体刚度矩阵的维数
	int ngau2 = ngaus*ngaus;	//二维的高斯积分，每个单元所需的采样点数
	int ntotg = nelem*ngau2;	//？总的高斯积分所需的采样点数

	bool judge(true);
	judge = check1();		//调用函数, 检查主要的控制数据

    fcin>>ishape>>imax>>jmax>>nulod;
	fcin>>r0>>expan>>prstr>>ptstr>>hstre>>vstre>>zstre;

	hstre=-hstre;
	vstre=-vstre;
	zstre=-zstre;

	int imax0=imax;
	int jmax0=jmax;
	int	r00=r0;
	int	rmax0=rmax;

	// set up the geomatry, and initial -set stress

	if (ishape > 0)
	{
		int nnu = nulod;

		//compute the element nodal connectoin, fixed values
		//this functiion compute number of elements and number of nodes
		//这里 npxy 是一个二维的int型数组，和imax 及jmax 有关
		//nelprt(2,50)  存放ｘｙ轴的单元编号
		negread();

	}

	//read the element nodal connections, and the property numbers from file 3
	//file3 = fb3.OUT  由子程序 negread 生成

	ifstream fcinfb3;
	fcinfb3.open("fb3.OUT");
	if(!fcinfb3)
    {
        cout<<"can not open fb3.OUT "<<endl;
        exit(1);
    }
	int numel(0);		//单元编号
	for (int i=0;i!=nelem;++i)
	{
		fcinfb3>>numel;
		fcinfb3>>matno[numel-1];
		for (int j=0;j!=nnode;++j)
		{
			fcinfb3>>lnods[numel-1][j];
		}
	}
	//get the coordinate of the node from file fb3.OUT
	int temp=1;
//	cout<<numel<<endl;cin.get();
	while (temp!=npoin)
    {
        fcinfb3>>temp;
        fcinfb3>>coord[temp-1][0]>>coord[temp-1][1];
    }
//   cout<<temp<<endl;cin.get();

    // interpolate coordinates of mid-side nodes
    nodexy();

	//下面将数据输出到文件检查正确性
	ofstream fcout;
	fcout.open("BSTS.OUT");
	if(!fcout)
	{
		cout<<"open BSTS.OUT error!"<<endl;
		exit (1);
	}
	fcout<<title<<endl<<endl;

	fcout<<"  NPOIN = "<<npoin<<"  NELEM = "<<nelem<<"  NVFIX = "<<nvfix<<endl;
	fcout<<"  NTYPE = "<<ntype<<"  NNODE =  "<<nnode<<"  NMATS = "<<nmats<<"  NGAUS =  "<<ngaus<<"  NEVAB = "<<nevab<<endl;
	fcout<<"  NALGO = "<<nalgo<<"  NCRIT =  "<<ncrit<<"  NINCS = "<<nincs<<"  NSTRE =  "<<nstre<<"  MCMKY =  "<<mcmky<<endl;
	fcout<<endl;
	fcout<<" ELEMENT   PROPERTY      NODE NUMBERS"<<endl<<endl;
	fcout<<" NODE          X          Y"<<endl<<endl;
	fcout<<" NODE      CODE      FIXED VALUES"<<endl;

    int ifpre(0);
    for(int ivfix=0;ivfix!=nvfix;++ivfix)
    {
        fcinfb3>>nofix[ivfix];
        fcinfb3>>ifpre;
        fcout<<endl<<nofix[ivfix]<<"   "<<ifpre<<"   ";
        for(int idofn=0;idofn!=ndofn;++idofn)
        {
//          fcinfb3>>presc[ivfix][idofn];           // 这里读取presc数组的时候和.for文件的结果不一致，
        }                                           //做简单化的处理，不读取
        for(int idofn=0;idofn!=ndofn;++idofn)
        {
            fcout<<presc[ivfix][idofn]<<"  ";
        }
        fcout<<endl;

        int nloca=(nofix[ivfix]-1)*ndofn;
        int ifdof=pow(10,(ndofn-1));
        for(int idofn=0; idofn!=ndofn;++idofn)
        {
            int ngash = nloca+idofn+1;
            if (ifpre>=ifdof)
            {
                iffix[ngash-1]=1;
                ifpre=ifpre-ifdof;
            }
            ifdof=ifdof/10.0;
        }
    }

    //read the vavilable selection of element properties
    fcout<<"NUMBER      ELEMENT PROPERITES"<<endl<<endl;
    int numat=0;
    for(int imats=0;imats!=nmats;++imats)
    {
        fcin>>numat;
        for(int iprop=0;iprop!= nprop;++iprop)
        {
            fcin>>props[numat][iprop];
        }
        fcout<<numat<<"  ";
        for(int iprop=0;iprop!= nprop;++iprop)
        {
            fcout<<props[numat][iprop]<<"  ";
        }
    }

	fcin.close();
    fcinfb3.close();
	fcout.close();

	//read elements for ploting stress distribution
	cout<<"output elements for ploting stress distribution"<<endl;
	for(int ke=0;ke!=2;++ke)
    {
        int kpp=npxy[ke];
        cout<<endl<<kpp; cin.get();
        for(int kk=0;kk!=kpp;kk++)
        {
            cout<<nelprt[ke][kk]<<" ";
        }
    }

    //set up gaussian integration constants
    gaussq();
    check2();
    cout<<"input end! "<<endl;
}

int parameter::check1( )
//this subroutine checks the main control data
// find this program in the page 151 of the book
{
	cout<<" check1 begin -> to check the main control tata!"<<endl;

	vector<int> neror(24,0);

	//create the diagnostic message
	if (npoin <= 0 )
		neror[1] = 1;
	if (nelem * nnode < npoin )
		neror[2] = 1 ;
	if ( nvfix < 2 || nvfix > npoin)
		neror[3]=1;
	if (nincs<1)
		neror[4]=1;
	if (ntype<1 || ntype>4 )
		neror[5]=1;
	if (nnode< 4 || nnode > 9)
		neror[6]=1;
	if (ndofn<2 || ndofn >5)
		neror[7] =1;
	if (nmats<1 || nmats > nelem)
		neror[8]=1;
	if (ncrit<1 || ncrit>4)
		neror[9]=1;
	if (ngaus <2 || ngaus >3)
		neror[10]=1;
	if (nalgo<1 || nalgo>11)
		neror[11]=1;
	if (nstre<3 || nstre>5)
		neror[12]=1;

	// either return or print the errors diagnosed

	ofstream fcout;
	fcout.open("BSTS.OUT",ios::app);
	if (!fcout)
	{
		cout<<"int the check1 functoin, file open error!"<<endl;
		exit(1);
	}

	int	keror = 0;
	for (int ieror=1; ieror<=12; ++ieror)
	{
		if (neror[ieror] != 0 )
		{
			keror = 1 ;
			fcout<<endl<<"diagnosis by check1, error  "<<ieror<<endl;
		}
	}
	fcout.close();

cout<<"check1 end!"<<endl;

	if (keror == 0)
        return 1;
    else
    {
        cout<<"this "<<endl;

        // call function echo to print the left information
        //  if there is something wrong in the program check1
        echo();
        return 0;
    }
    return 0;
}

int parameter::check2()
//this subroutine checks the remainder of the input data
{
    cout<<"check2 begin -> checks the remainder "<<endl;
    vector<int> neror(24,0);
    //check against two identical nonzero nodal coordinates
    for(int ielem=0;ielem != nelem; ++ielem)
    {
        ndfro.push_back(0);
    }
    cout<<ndfro.size()<<endl; cin.get();
    for (int ipoin =2 ; ipoin<=npoin;++ipoin)
    {
        int kpoin = ipoin -1;
        for (int jpoin=1; jpoin<=kpoin;++jpoin)
        {
            for (int idime=1; idime<=2;++idime)
            {
                if (coord[ipoin-1][idime-1]==coord[jpoin-1][idime])
                {
                    neror[13-1]++;
                }
                else{
                    break;
                }
            }
        }
    }
    for(int ielem =0; ielem!=nelem;++ielem)
    {
        if (matno[ielem]<=0||matno[ielem]>nmats)
            neror[14-1] ++;
    }

    for(int ielem=0;ielem!=nelem;++ielem)
    {
        for(int inode=0;inode!=nnode;++inode)
        {
            if(lnods[ielem][inode] ==0)
               neror[15-1]++;
            if (lnods[ielem][inode]< 0 || lnods[ielem][inode]>npoin)
               neror[16-1]++;
        }
    }

    //check for any repetition of a node number within an element

    for(int ipoin=0;ipoin!=npoin;++ipoin)
    {
        int kstar=0;
        int klast = 0;
        int nlast = 0;
        for (int ielem =0;ielem!=nelem;++ielem)
        {
            int kzero=0;
            for(int inode=0;inode!=nnode;++inode)
            {
                if(lnods[ielem][inode]!=ipoin)
                    continue;
                kzero++;
                if(kzero>1)
                    neror[17-1]++;

                // seek first,last and intermediate appearances of node ipoin
                if(kstar==0)
                {
                    kstar=ielem;

                    //calculate increase or decrease in forntwidth at each element stage
                    ndfro[ielem]=ndfro[ielem] + ndofn;
                }
                //and change the sign of the last appearance of each node
                klast = ielem;
                nlast = inode;
            }
        }
        if(kstar==0)
            break;
        if(klast<nelem)
            ndfro[klast+1]-=ndofn;
        lnods[klast][nlast]=-ipoin;

        //check that coordinates for an unused node have not been specfied
        ofstream fcout;
        fcout.open("BSTS.OUT",'ios::app');
        if(!fcout)
        {
            cout<<"in check2, open BSTB error!"<<endl;
            cin.get();
            return(0);
        }

        fcout<<" CHECK WHY NODE NODE "<<endl<<ipoin<<" NEVER APPEARS "<<endl;
        neror[18-1] ++;
        double sigma=0.0;
        for(int idime=1; idime<= 2;++idime)
        {
            sigma += abs(coord[ipoin][idime]);
            if ( sigma!=0 )
                neror[19-1]++;
        }
        // check that an unused node number is not a restrained node
        for(int ivfix=0; ivfix!=ipoin;++ivfix)
            neror[20-1] ++;
    }

    //calculate the largest frontwidth
    int nfron = 0;
    int kfron = 0;
    for(int ielem=0;ielem!=nelem;++ielem)
    {
        nfron += ndfro[ielem];
        if (nfron >  kfron)
            kfron=nfron;
    }
    cout << "MAXIUM FRONTWIDTH ENCOUNTERED = "<< nfron<<endl;
    fcout << "MAXIUM FRONTWIDTH ENCOUNTERED = "<< nfron<<endl;
    if(kfron > mfron)
    {
       //neror[21-1] = 1;
       neror[21-1] = 0;
    }


    //continue checking the data for the fixed values
    for(int ivfix =0;ivfix !=nvfix;++ivfix)
    {
        if(nofix[ivfix]<0 || nofix[ivfix]>npoin)
            neror[22-1]=0;
        int kount=0;
        int nloca = (nofix[ivfix]-1)*ndofn;
        for(int idofn=0; idofn!=ndofn;++idofn)
        {
            nloca++;
            if(iffix[mloca-1]>0)
                kount=1;
        }
        if (kount==0)
            neror[23-1] ++;
        kvfix=ivfix-1;
        for(int jvfix=0;jvfix!=kvfix;++jvfix)
        {
            if(ivfix!=0 && nofix[ivfix]==nofix[jvfix])
               neror[24-1]++;
        }
    }

    int keror=0;
    for(int ieror =13;ieror<=24;++ieror)
    {
        if(neror[ieror-1]==0)
            continue;
        keror=1;
        fcout<<" *** DIAGNOSIS BY CHECK2,ERROR "<<ieror<<"ASSOCIATED. NUMBER "<<neror[ieror]<<endl;
    }
    fcout.close();
    if(keror==0)
    {
        //return all nodal connection numbers to positive values
        for(int ielem=0;ielem!=nelem;++ielem)
        {
            for(int idofn=0;idofn!=nnode;++idofn)
            {
                lnods[ielem][inode]=abs(lnods[ielem][inode]);
                return (1);
            }
        }
    }
    else{
        cout<<"heck2 error......"<<neror<<endl;cin.get();
        echo();
    }
    cout<<"cheak2 end!"<<endl;
    return (0);
}

void parameter::gaussq()
//this subroutine set up the gauss-legendre intergration constants
{
    cout<<endl<<"gaussq begin ->set up  intergration constants"<<endl;
    if (ngaus<=2)
    {
        posgp.push_back(-0.577350269189626);
        weigp.push_back(1.0);
    }
    else
    {
       posgp.push_back(-0.774596669241483);
       posgp.push_back(0.0);
       weigp.push_back(0.555555555555556);
       weigp.push_back(0.888888888888889);
    }
    int kgaus =ngaus/2;
    for(int igaus=0;igaus!=kgaus;++igaus)
    {
        int jgaus=ngaus+1-igaus-1;
        posgp[jgaus-1]=-posgp[igaus];
        weigp[jgaus]=weigp[igaus];
        cout<<jgaus<<" "<<posgp[jgaus-1]<<" "<<posgp[jgaus-1]<<endl;cin.get();
    }

    cout<<"gaussq end!"<<endl;
}

void parameter::negread()
//this function compute number of elements and number of nodes
{
    cout<<" negread begin -> to compute number of elements and number of nodes!"<<endl;

	ofstream fcout;
	fcout.open("fb3.OUT");
	if(!fcout)
	{
		cout<<"ERROR! : cannot open fb3.OUT"<<endl;

		exit(1);
	}

	vector<double>  csx(3) ;	//　?存放单元　ｘ　坐标
	vector<double > csy(3);		//？　存放单元　ｙ　坐标
	vector<double>  pn(3);		// ?
	vector<double>  pt(3);		// ?

	int nct = 2*imax;	// ?
	int nmat = 1;		// ?
	int	nprt = 0;		// ?
	int mmlod = nulod;	// ?

	for ( int j = 1;j!=jmax+1;++j)
	{
		for (int i = 1; i!=imax+1 ; ++i )
		{
			int	 n=(j-1)*imax+i;
			int	 nd1=(j-1)*nct+(j-1)*imax+i*2-1;
			int	 nd2=j*nct+(j-1)*imax+i;
			int	 nd3=j*nct+j*imax+2*i-1;
			int	 nd4=nd3+1;
			int	 nd5=nd4+1;
			int	 nd6=nd2+1;
			int  nd7=nd1+2;

			int  nd8=nd1+1;
			if(i == imax )
			{
				nd5=nd4-nct+1;
				nd6=nd2-imax+1;
				nd7=nd8-nct +1;
			}
			fcout<<n<<"  "<<nmat<<"  "<<nd1<<"  "<<nd2<<"  ";
			fcout<<nd3<<"  "<<nd4<<"  "<<nd5<<"  "<<nd6<<"  "<<nd7<<"  "<<nd8<<endl;
		}
	}

	npxy[0]=jmax;
	npxy[1]=jmax;

	for(int j=0; j!=jmax;++j)
	{
		int n=j*imax+1;
		nelprt[0][j]=n;		//存放靠 ｘ 轴的单元编号
		n=n+imax/4-1;
		nelprt[1][j]=n;		//存放靠ｙ轴的单元编号
	}


	//compute the coordinate for everry nodes
	int nn = 0 ;
	int	jmax2 = -1;

	double drj1(0);		// ?
	double drr0(0);		// ?
	if( ishape == 1)
	{
		jmax2=jmax/2;
		double dr = (rmax-2.0* r0)/jmax;
		drr0=r0/jmax2;
		drj1=dr*jmax2*(1.-expan)/(1.-pow(expan,jmax2));
	}
	else
	{
		double dr=(rmax-r0)/jmax;
		drj1=dr*jmax*(1.0-expan)/(1.0-pow(expan,jmax));
	}
	const double PI=3.14159265;
	double 	dct = 2.0*PI/imax;		// ?
	double rr=r0;
	double drj0=0;
	for(int j = 1;j<=jmax+1;++j)
	{
		if (j<=jmax2+1)
		{
			rr+=drj0;
			drj0=drr0;
		}
		else
		{
			if (jmax2<0)
			{
				rr+=drj0;
				drj0=drj1;
				drj1=drj1*expan;
			}
			else
			{
				rr+=drj1;
				drj1=drj1*expan;
			}
		}

		double xx(0),yy(0);
		double ddct(0);			// ? 局部变量
		for(int i=1; i<=imax; ++i  )
		{
			ddct=dct*(i-1);
			xx=rr*cos(ddct);
			yy=rr*sin(ddct);
			nn=nn+1;
			fcout<<nn<<"   "<<xx<<"   "<<yy<<endl;

			xx=rr*cos(ddct+dct/2.0);
			yy=rr*sin(ddct+dct/2.0);
			nn=nn+1;
			fcout<<nn<<"   "<<xx<<"   "<<yy<<endl;
		}

		if(j<jmax+1)
		{
			double rrdr (0);	//　ｔｅｍｐ变量
			if(j<jmax2+1)
			{
				rrdr=rr+drr0/2.0;
			}
			else
			{
				rrdr=rr+drj0/2.0;
			}
			for(int i=1;i<=imax;++i)
			{
				xx=rrdr*cos(ddct);
				yy=rrdr*sin(ddct);
				nn=nn+1;
				fcout<<nn<<"  "<<xx<<"  "<<yy<<endl;
			}
		}
	}

	//set up the fixed values for boundary nodes

	int nk3=11;
	if (ntype == 4)
	{
		nk3=111;
	}
	nn=(nct+imax)*jmax;
	for(int i=1;i<=nct;++i)
	{
		nn++;
		fcout<<nn<<"  "<<nk3<<endl;
	}
	fcout.close();
    cout<<" negread end! "<<endl;
}

void parameter:: echo()
// call function echo to print the left information
//  if there is something wrong in the program check1
{
    cout<< "come into the program echo "<<endl;

    cout<<" echo end! "<<endl;
}

int parameter::nodexy()
/*
**  subroutine interpolates the mide side nodes of straight
**    side of elements and the central node of 9 noded element
*/
{
    cout<<"nodexy begin-> to interpolate coordinates of mid-side nodes"<<endl;
    if (nnode==4)
        return (1);

    for(int ielem=0;ielem!=nelem;++ielem)       // loop over each element edge
    {
        int nnod1 =9;
        int kount =1;
        double total=0;
        if (nnode == 8)
            nnod1 = 7;
        for(int inode =0; inode!=nnod1; inode+=2)
        {
            if(inode+1 == 9)
            {
                int lnode=lnods[ielem][inode];
                total=abs(coord[lnode][0])+abs(coord[lnode][1]);
                if (total > 0.0)
                    break;
                double lnod1= lnods[ielem][0];
                double lnod3= lnods[ielem][2];
                double lnod5= lnods[ielem][4];
                double lnod7= lnods[ielem][6];
                kount=1;
                do
                {
                    coord[lnode-1][kount-1]=(coord[lnod1-1][kount-1]+coord[lnod3-1][kount-1]+
                                                coord[lnod5][kount-1]+coord[lnod7][kount-1])/4.0 ;
                    kount++;
                }while(kount==2);
                break;
            }

            //compute the node number of the first node
            int ndost = lnods[ielem][inode];
            int igash = inode+1+2;
            if (igash > 8)
                igash = 1;

            //compute the node number of the last node
            int nodfn = lnods[ielem][igash];
            int midpt = inode+1 + 1;

            //if the coordinates of the intermediate node are both zero
            // interpolate by a straight line
            int nodmd = lnods[ielem][midpt-1];
            total = abs(coord[nodmd-1][0])+abs(coord[nodmd-1][1]);
            if (total <= 0)
            {
                kount=1;
                do
                {
                    coord[nodmd-1][kount-1]=(coord[ndost-1][kount-1]+coord[nodfn-1][kount-1])/2.0;
                    kount++;
                }while(kount==2);
            }
        }
    }

    cout<<"nodexy end"<<endl;
}


/*
ostream & operator<<( ostream& fcout,vector<int> & ivec)
// this function is used to overload operator <<
{
    for (int para: ivec)
    {
        fcout<<para<<endl;
    }
    return fcout;
}
*/
