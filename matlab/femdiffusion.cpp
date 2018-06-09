/*
 * =============================================================
 * FEM Kernel for diffusioin equation
 * Author: Qianqian Fang(fangq < at > nmr.mgh.harvard.edu)
 * Date: 2007/11/21
 * Version: 0.5.0
 *
 * Model:
 *    [A]*c=-OC
 * where [A]=<D delPhi_i,delPhi_j>
 * OC is the vector for oxygen consumption at each node
 * boundary condition is dot(delc,normal)=0
 * =============================================================
 */

#include <string.h>
#include <exception>
#include <stdio.h>

#include "mex.h"
#include "femdiffusion.h"

#define N137_REFF 0.468407140351385

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
  Config cfg;
  tetmesh mesh;
  Forward femdata;
  Jacobian jac;
  int isjacobian=0;
  int nfields, ifield;

  if (nrhs<1){
     rb3_usage();
     return;
  }

  /* Get the matrices from the input list */
  
  if (!mxIsStruct(prhs[0]))
     mexErrMsgTxt("Input must be a structure.");

  try{
    config_init(&cfg);
    mesh_init(&mesh);
    forward_init(&femdata,&mesh);
    jacobian_init(&jac);
    
    nfields = mxGetNumberOfFields(prhs[0]);
    for (ifield = 0; ifield < nfields; ifield++) { /* how many input struct fields */
      mxArray *tmp = mxGetFieldByNumber(prhs[0], 0, ifield);
      if (tmp == NULL) 
            continue;
      mcx_set_field(prhs[0],tmp,ifield,&cfg,&mesh,&jac);
    }

    if(mesh.nn<=0||mesh.ne<=0)
      mexErrMsgTxt("dimensions of the input or output parameters are incorrect.");

    if(nrhs>1){
        const dimtype *arraydim;
	arraydim=mxGetDimensions(prhs[1]);
	if(arraydim[1]<3 || arraydim[0]==0)
	    MEXERROR("sd must have at least 3 columns");
        jac.sd=mxGetPr(prhs[1]);
	jac.nsd=arraydim[0];
	jac.nsdcol=arraydim[1];
	isjacobian=nlhs;
        if(nrhs<=2)
	    MEXERROR("sd and phi must be given to calculate Jacobians");
	arraydim=mxGetDimensions(prhs[2]);
	if(arraydim[0]!=mesh.nn)
	    MEXERROR("phi's row number must be the same as the number of node");
	jac.Phir=mxGetPr(prhs[2]); // phi: [mesh.nn x noptodes]
	if(mxIsComplex(prhs[2]))
	    jac.Phii=mxGetPi(prhs[2]); // phi: [mesh.nn x noptodes]

        if(nrhs>3){
	    arraydim=mxGetDimensions(prhs[3]);
	    if(arraydim[0]!=mesh.ne || arraydim[1]!=10)
	        MEXERROR("deldotdel must have be a 2D array with size: elem number x 10");
	    jac.deldotdel=mxGetPr(prhs[3]); // phi: [mesh.nn x jac.nsd]
	}
	if(nrhs>4){
	    jac.isnodal=((int)*(mxGetPr(prhs[4])));
	}
    }

    /* Assign a pointer to the output. */  
    if(cfg.omega>0.f){
        if(nlhs>=1){
	    if(isjacobian==0){
                plhs[0] = mxCreateDoubleMatrix(1,mesh.nn, mxCOMPLEX);
                femdata.Di= mxGetPi(plhs[0]);
	    }else{
                plhs[0] = mxCreateDoubleMatrix(jac.nsd, (jac.isnodal ? mesh.nn: mesh.ne), mxCOMPLEX);
                jac.Jmuai = mxGetPi(plhs[0]);
	    }
	}
	if(nlhs>=2){
	    if(isjacobian==0){
                plhs[1] = mxCreateDoubleMatrix(1,mesh.ntot, mxCOMPLEX);
                femdata.Ai= mxGetPi(plhs[1]);
	    }else{
                plhs[1] = mxCreateDoubleMatrix(jac.nsd, (jac.isnodal ? mesh.nn: mesh.ne), mxCOMPLEX);
                jac.Jdi = mxGetPi(plhs[1]);
	    }
	}
    }else{
        if(isjacobian==0){
          if(nlhs>=1)
            plhs[0] = mxCreateDoubleMatrix(1,mesh.nn, mxREAL);
	  if(nlhs>=2)
            plhs[1] = mxCreateDoubleMatrix(1,mesh.ntot, mxREAL);
	}else{
          if(nlhs>=1)
            plhs[0] = mxCreateDoubleMatrix(jac.nsd, (jac.isnodal ? mesh.nn: mesh.ne), mxREAL);
	  if(nlhs>=2)
            plhs[1] = mxCreateDoubleMatrix(jac.nsd, (jac.isnodal ? mesh.nn: mesh.ne), mxREAL);
	}
    }
    if(isjacobian==0){
      if(nlhs>=1)
        femdata.Dr= mxGetPr(plhs[0]);
      if(nlhs>=2)
        femdata.Ar= mxGetPr(plhs[1]);
    }else{
      if(nlhs>=1)
        jac.Jmuar = mxGetPr(plhs[0]);
      if(nlhs>=2)
        jac.Jdr = mxGetPr(plhs[1]);
    }
    if(nlhs>=3){
	plhs[2] = mxCreateDoubleMatrix(10,mesh.ne, mxREAL);
        femdata.deldotdel = mxGetPr(plhs[2]);
    }
    
    if(isjacobian==0){
	femdiffusion(&cfg,&mesh,&femdata);
	femdiffusion_boundary(&cfg,&mesh,&femdata);
    }else{
        femjacobian(&cfg,&mesh,&jac);
    }

    /** \subsection sclean End the simulation */
    mesh_clear(&mesh);
    //mcx_clearcfg(&cfg);

  }catch(const char *err){
      mexPrintf("Error: %s\n",err);
  }catch(const std::exception &err){
      mexPrintf("C++ Error: %s\n",err.what());
  }catch(...){
      mexPrintf("Unknown Exception");
  }
}

void mcx_set_field(const mxArray *root,const mxArray *item,int idx, Config *cfg, tetmesh *mesh, Jacobian *jac){
    const char *name=mxGetFieldNameByNumber(root,idx);
    const dimtype *arraydim;
    int i,j;

    GET_1ST_FIELD(cfg,omega)
    GET_VEC34_FIELD(cfg,srcpos)
    GET_VEC34_FIELD(cfg,srcdir)
    GET_ONE_FIELD(mesh,e0)
    GET_ONE_FIELD(mesh,isreoriented)
    GET_ONE_FIELD(cfg,reff)
    else if(strcmp(name,"node")==0){
        arraydim=mxGetDimensions(item);
	if(arraydim[0]<=0 || arraydim[1]!=3)
            MEXERROR("the 'node' field must have 3 columns (x,y,z)");
        double *val=mxGetPr(item);
        mesh->nn=arraydim[0];
	if(mesh->node) free(mesh->node);
        mesh->node=(float3 *)calloc(sizeof(float3),mesh->nn);
        for(j=0;j<3;j++)
          for(i=0;i<mesh->nn;i++)
             ((double *)(&mesh->node[i]))[j]=val[j*mesh->nn+i];
        printf("rbm.nn=%d;\n",mesh->nn);
    }else if(strcmp(name,"elem")==0){
        arraydim=mxGetDimensions(item);
	if(arraydim[0]<=0 || arraydim[1]!=4)
            MEXERROR("the 'elem' field must have 4 columns (e1,e2,e3,e4)");
        double *val=mxGetPr(item);
        mesh->ne=arraydim[0];
	if(mesh->elem) free(mesh->elem);
        mesh->elem=(int4 *)calloc(sizeof(int4),mesh->ne);
        for(j=0;j<4;j++)
          for(i=0;i<mesh->ne;i++)
             ((int *)(&mesh->elem[i]))[j]=(int)val[j*mesh->ne+i]-1;
        printf("rbm.ne=%d;\n",mesh->ne);
    }else if(strcmp(name,"face")==0){
        arraydim=mxGetDimensions(item);
	if(arraydim[0]<=0 || arraydim[1]!=3)
            MEXERROR("the 'face' field must have 3 columns (e1,e2,e3)");
        double *val=mxGetPr(item);
        mesh->nf=arraydim[0];
	if(mesh->face) free(mesh->face);
        mesh->face=(int3 *)calloc(sizeof(int3),mesh->nf);
        for(j=0;j<3;j++)
          for(i=0;i<mesh->nf;i++)
             ((int *)(&mesh->face[i]))[j]=(int)val[j*mesh->nf+i]-1;
        printf("rbm.nf=%d;\n",mesh->nf);
    }else if(strcmp(name,"elemprop")==0){
        arraydim=mxGetDimensions(item);
	if(MAX(arraydim[0],arraydim[1])==0)
            MEXERROR("the 'elemprop' field can not be empty");
        double *val=mxGetPr(item);
        mesh->ne=MAX(arraydim[0],arraydim[1]);
	if(mesh->type) free(mesh->type);
	mesh->type=(int  *)malloc(sizeof(int )*mesh->ne);
        for(i=0;i<mesh->ne;i++)
           mesh->type[i]=val[i];
        printf("rbm.ne=%d;\n",mesh->ne);
    }else if(strcmp(name,"evol")==0){
        arraydim=mxGetDimensions(item);
	if(MAX(arraydim[0],arraydim[1])==0)
            MEXERROR("the 'evol' field can not be empty");
        double *val=mxGetPr(item);
        mesh->ne=MAX(arraydim[0],arraydim[1]);
	if(mesh->evol) free(mesh->evol);
        mesh->evol=(double *)malloc(sizeof(double)*mesh->ne);
        for(i=0;i<mesh->ne;i++)
           mesh->evol[i]=val[i];
        printf("rbm.evol=%d;\n",mesh->ne);
    }else if(strcmp(name,"area")==0){
        arraydim=mxGetDimensions(item);
	if(MAX(arraydim[0],arraydim[1])==0)
            MEXERROR("the 'evol' field can not be empty");
        double *val=mxGetPr(item);
        mesh->nf=MAX(arraydim[0],arraydim[1]);
	if(mesh->area) free(mesh->area);
        mesh->area=(double *)malloc(sizeof(double)*mesh->nf);
        for(i=0;i<mesh->nf;i++)
           mesh->area[i]=val[i];
        printf("rbm.area=%d;\n",mesh->nf);
    }else if(strcmp(name,"idxcount")==0){
        arraydim=mxGetDimensions(item);
	if(MAX(arraydim[0],arraydim[1])==0)
            MEXERROR("the 'evol' field can not be empty");
        double *val=mxGetPr(item);
        mesh->nn=MAX(arraydim[0],arraydim[1]);
	if(mesh->idxcount) free(mesh->idxcount);
        mesh->idxcount=(int *)malloc(sizeof(int)*mesh->nn);
        for(i=0;i<mesh->nn;i++)
           mesh->idxcount[i]=val[i];
        printf("rbm.idxcount=%d;\n",mesh->nn);
    }else if(strcmp(name,"idxsum")==0){
        arraydim=mxGetDimensions(item);
	if(MAX(arraydim[0],arraydim[1])==0)
            MEXERROR("the 'evol' field can not be empty");
        double *val=mxGetPr(item);
        mesh->nn=MAX(arraydim[0],arraydim[1]);
	if(mesh->idxsum) free(mesh->idxsum);
        mesh->idxsum=(int *)malloc(sizeof(int)*mesh->nn);
        for(i=0;i<mesh->nn;i++)
           mesh->idxsum[i]=val[i];
        printf("rbm.idxsum=%d;\n",mesh->nn);
    }else if(strcmp(name,"rows")==0){
        arraydim=mxGetDimensions(item);
	if(MAX(arraydim[0],arraydim[1])==0)
            MEXERROR("the 'rows' field can not be empty");
        double *val=mxGetPr(item);
        mesh->ntot=MAX(arraydim[0],arraydim[1]);
	if(mesh->rows) free(mesh->rows);
        mesh->rows=(int *)malloc(sizeof(int)*mesh->ntot);
        for(i=0;i<mesh->ntot;i++)
           mesh->rows[i]=(int)val[i]-1;
        printf("rbm.rows=%d;\n",mesh->ntot);
    }else if(strcmp(name,"cols")==0){
        arraydim=mxGetDimensions(item);
	if(MAX(arraydim[0],arraydim[1])==0)
            MEXERROR("the 'cols' field can not be empty");
        double *val=mxGetPr(item);
        mesh->ntot=MAX(arraydim[0],arraydim[1]);
	if(mesh->cols) free(mesh->cols);
        mesh->cols=(int *)malloc(sizeof(int)*mesh->ntot);
        for(i=0;i<mesh->ntot;i++)
           mesh->cols[i]=(int)val[i]-1;
        printf("rbm.cols=%d;\n",mesh->ntot);
    }else if(strcmp(name,"prop")==0){
        arraydim=mxGetDimensions(item);
        if(arraydim[0]>0 && arraydim[1]!=4)
            MEXERROR("the 'prop' field must have 4 columns (mua,mus,g,n)");
        double *val=mxGetPr(item);
        mesh->prop=arraydim[0]-1;
        if(mesh->med) free(mesh->med);
        mesh->med=(medium *)calloc(sizeof(medium),mesh->prop+1);
        for(j=0;j<4;j++)
          for(i=0;i<=mesh->prop;i++)
             ((double *)(&mesh->med[i]))[j]=val[j*(mesh->prop+1)+i];
        printf("rbm.prop=%d;\n",mesh->prop);
    }else if(strcmp(name,"deldotdel")==0){
        if(jac->deldotdel==NULL){
	    arraydim=mxGetDimensions(item);
	    if(arraydim[1]!=10)
	        MEXERROR("deldotphi must have be a 2D array with size: node number x 10");
	    jac->deldotdel=mxGetPr(item); // phi: [mesh.nn x jac.nsd]
            printf("rbm.deldotdel=<%d,10>;\n",arraydim[0]);
	}
    }else if(strcmp(name,"detpos")==0){
        // do nothing
    }else if(strcmp(name,"nvol")==0){
        // do nothing
    }else{
        printf("WARNING: redundant field '%s'\n",name);
    }
    mexEvalString("pause(.0001);");
}

void config_init(Config *cfg){
	cfg->omega=0.f;
	cfg->reff=N137_REFF;
	memset(&(cfg->srcpos),0,sizeof(float4));
	memset(&(cfg->srcdir),0,sizeof(float4));
}

void mcx_clearcfg(Config *cfg){
}

void forward_init(Forward *fem, tetmesh *m){
	fem->mesh=m;
	fem->Ar=NULL;
	fem->Ai=NULL;
	fem->Dr=NULL;
	fem->Di=NULL;
	fem->deldotdel=NULL;
}

void jacobian_init(Jacobian *jac){
	jac->isnodal=1;
	jac->nsd=0;
	jac->nsdcol=0;
	jac->Jmuar=NULL;
	jac->Jmuai=NULL;
	jac->Jdr=NULL;
	jac->Jdi=NULL;
	jac->Phir=NULL;
	jac->Phii=NULL;
	jac->sd=NULL;
	jac->deldotdel=NULL;
}

void mesh_init(tetmesh *mesh){
	mesh->nn=0;
	mesh->ne=0;
	mesh->nf=0;
	mesh->prop=0;
	mesh->node=NULL;
	mesh->elem=NULL;
	mesh->face=NULL;
	mesh->type=NULL;
	mesh->med=NULL;
	mesh->evol=NULL;
	mesh->area=NULL;
	mesh->rows=NULL;
	mesh->cols=NULL;
	mesh->idxcount=NULL;
	mesh->idxsum=NULL;
}

void mesh_clear(tetmesh *mesh){
	mesh->nn=0;
	mesh->ne=0;
	mesh->nf=0;
        mesh->ntot=0;
        mesh->prop=0;
	if(mesh->node){
		free(mesh->node);
		mesh->node=NULL;
	}
	if(mesh->elem){
		free(mesh->elem);
		mesh->elem=NULL;
	}
	if(mesh->face){
		free(mesh->face);
		mesh->face=NULL;
	}
	if(mesh->area){
		free(mesh->area);
		mesh->area=NULL;
	}
	if(mesh->type){
		free(mesh->type);
		mesh->type=NULL;
	}
	if(mesh->med){
		free(mesh->med);
		mesh->med=NULL;
	}
        if(mesh->evol){
                free(mesh->evol);
                mesh->evol=NULL;
        }
        if(mesh->rows){
                free(mesh->rows);
                mesh->rows=NULL;
        }
        if(mesh->cols){
                free(mesh->cols);
                mesh->cols=NULL;
        }
        if(mesh->idxcount){
                free(mesh->idxcount);
                mesh->idxcount=NULL;
        }
        if(mesh->idxsum){
                free(mesh->idxsum);
                mesh->idxsum=NULL;
        }
}

void femjacobian(Config *cfg,tetmesh *mesh, Jacobian *jac){
    int t, sd, sid, rid, i;
    int pairs[3][10]={{0,0,0,1,1,2},{1,2,3,1,2,3},{1,2,3,5,6,8}}; //pairs[0]<->[1], local node pairs, [3] pos in deldotdel
    int inode[4]={0,4,7,9}; // position in deldotdel for diagonal

    if(jac->Jmuar)
       memset(jac->Jmuar,0,jac->nsd*(jac->isnodal?mesh->nn:mesh->ne)*sizeof(double));
    if(jac->Jmuai)
       memset(jac->Jmuai,0,jac->nsd*(jac->isnodal?mesh->nn:mesh->ne)*sizeof(double));
    if(jac->Jdr)
       memset(jac->Jdr,0,jac->nsd*(jac->isnodal?mesh->nn:mesh->ne)*sizeof(double));
    if(jac->Jdi)
       memset(jac->Jdi,0,jac->nsd*(jac->isnodal?mesh->nn:mesh->ne)*sizeof(double));

    // loop over all elements
    for(t=0;t<mesh->ne;t++){
        int *ee=(int *)(mesh->elem+t);

        // first calculate real-part of Jmua and Jd
        for(sd=0;sd<jac->nsd;sd++){  // loop over sd pairs
	    double jdr=0.0, jmuar=0.0, dtemp;
            sid=(int)jac->sd[sd]-1;
	    rid=(int)jac->sd[sd+jac->nsd]-1;

	    // accummulate diagonal terms
	    for(i=0;i<4;i++){
	        dtemp=jac->Phir[sid*mesh->nn+ee[i]]*jac->Phir[rid*mesh->nn+ee[i]];
		if(jac->Phii)
		    dtemp-=jac->Phii[sid*mesh->nn+ee[i]]*jac->Phii[rid*mesh->nn+ee[i]];
		jdr+=jac->deldotdel[inode[i]*mesh->ne+t]*dtemp;
		jmuar+=dtemp;
	    }
	    // accummulate off-diagonal terms
	    for(i=0;i<6;i++){
	        dtemp=jac->Phir[sid*mesh->nn+ee[pairs[0][i]]]*jac->Phir[rid*mesh->nn+ee[pairs[1][i]]] +
		      jac->Phir[rid*mesh->nn+ee[pairs[0][i]]]*jac->Phir[sid*mesh->nn+ee[pairs[1][i]]];
		if(jac->Phii)
		    dtemp-=jac->Phii[sid*mesh->nn+ee[pairs[0][i]]]*jac->Phii[rid*mesh->nn+ee[pairs[1][i]]] +
		           jac->Phii[rid*mesh->nn+ee[pairs[0][i]]]*jac->Phii[sid*mesh->nn+ee[pairs[1][i]]];
	        jmuar+=dtemp*0.5;
		jdr+=jac->deldotdel[pairs[2][i]*mesh->ne+t]*dtemp;
	    }
	    jmuar*=-0.1;
	    jdr=-jdr;
	    if(jac->isnodal){
	        for(i=0;i<4;i++)
		    jac->Jmuar[ee[i]*jac->nsd+sd]+=jmuar*0.25*mesh->evol[t];
	    }else{ // element-based jacobian
		jac->Jmuar[t*jac->nsd+sd]+=jmuar*mesh->evol[t];
	    }
	    if(jac->Jdr){
	        if(jac->isnodal){
	            for(i=0;i<4;i++)
	                jac->Jdr[ee[i]*jac->nsd+sd]+=jdr*0.25;
		}else{ // element-based jacobian
		    jac->Jdr[t*jac->nsd+sd]+=jdr;
	        }
	    }
	}
    }

    // if measurement has complex val, then calculate the imag part of J
    if(jac->Phii){
      for(t=0;t<mesh->ne;t++){
        int *ee=(int *)(mesh->elem+t);
    
        // first calculate real-part of Jmua and Jd
        for(sd=0;sd<jac->nsd;sd++){  // loop over sd pairs
	    double jdi=0.0, jmuai=0.0, dtemp;
            sid=(int)jac->sd[sd]-1;
	    rid=(int)jac->sd[sd+jac->nsd]-1;

	    // accummulate diagonal terms
	    for(i=0;i<4;i++){
	        dtemp=jac->Phir[sid*mesh->nn+ee[i]]*jac->Phii[rid*mesh->nn+ee[i]];
		if(jac->Phii)
		    dtemp+=jac->Phii[sid*mesh->nn+ee[i]]*jac->Phir[rid*mesh->nn+ee[i]];
		jdi+=jac->deldotdel[inode[i]*mesh->ne+t]*dtemp;
		jmuai+=dtemp;
	    }
	    // accummulate off-diagonal terms
	    for(i=0;i<6;i++){
	        dtemp=jac->Phir[sid*mesh->nn+ee[pairs[0][i]]]*jac->Phii[rid*mesh->nn+ee[pairs[1][i]]] +
		      jac->Phir[rid*mesh->nn+ee[pairs[0][i]]]*jac->Phii[sid*mesh->nn+ee[pairs[1][i]]];
		if(jac->Phii)
		    dtemp+=jac->Phii[sid*mesh->nn+ee[pairs[0][i]]]*jac->Phir[rid*mesh->nn+ee[pairs[1][i]]] +
		           jac->Phii[rid*mesh->nn+ee[pairs[0][i]]]*jac->Phir[sid*mesh->nn+ee[pairs[1][i]]];
	        jmuai+=dtemp*0.5;
		jdi+=jac->deldotdel[pairs[2][i]*mesh->ne+t]*dtemp;
	    }
	    jmuai*=-0.1;
	    jdi=-jdi;
	    if(jac->isnodal){ // nodal-based Jacobian
	        for(i=0;i<4;i++)
		    jac->Jmuai[ee[i]*jac->nsd+sd]+=jmuai*0.25*mesh->evol[t];
	    }else{ // element-based jacobian
		jac->Jmuai[t*jac->nsd+sd]+=jmuai*mesh->evol[t];
	    }
	    if(jac->Jdi){
	        if(jac->isnodal){
	            for(i=0;i<4;i++)
	                jac->Jdi[ee[i]*jac->nsd+sd]+=jdi*0.25;
		}else{ // element-based jacobian
		    jac->Jdi[t*jac->nsd+sd]+=jdi;
	        }
	    }
	}
      }
    }
}

void femdiffusion_boundary(Config *cfg,tetmesh *mesh, Forward *fem){
    int t;
    double val;
    double Reff=(1.0-cfg->reff)/(1.0+cfg->reff);

    for(t=0;t<mesh->nf;t++){
        val=mesh->area[t]*Reff*(1.0/12.0);
	fem->Dr[mesh->face[t].x]+=val;
	fem->Dr[mesh->face[t].y]+=val;
	fem->Dr[mesh->face[t].z]+=val;

	val=val*0.5;
	fem->Ar[sub2ind(mesh,mesh->face[t].x,mesh->face[t].y)]+=val;
	fem->Ar[sub2ind(mesh,mesh->face[t].x,mesh->face[t].z)]+=val;
	fem->Ar[sub2ind(mesh,mesh->face[t].y,mesh->face[t].z)]+=val;
    }
}

void femdiffusion(Config *cfg,tetmesh *mesh, Forward *fem){
    double mua, musp,dcoeff,nref;
    float3 *N;
    int i1,i2,i3,i4,i,j,t,ii,jj,ij;
    double derx[4],dery[4],derz[4];
    double ra,ri,Ve,RVe6,deldotdel,sm;
    
    /*initialize output variables*/
    if(fem->Dr)
	memset(fem->Dr,0,mesh->nn*sizeof(double));
    if(fem->Di)
	memset(fem->Di,0,mesh->nn*sizeof(double));
    if(fem->Ar)
	memset(fem->Ar,0,mesh->ntot*sizeof(double));
    if(fem->Ai)
	memset(fem->Ai,0,mesh->ntot*sizeof(double));

    /* map coordinate pointers*/
    N=mesh->node;
    
    /* loop over all elements*/
    for(t=0;t<mesh->ne;t++){
        RVe6=1.0/(mesh->evol[t]*6.0);
        Ve =mesh->evol[t];
        
        i1=mesh->elem[t].x;
        i2=mesh->elem[t].y;
        i3=mesh->elem[t].z;
        i4=mesh->elem[t].w;

	mua=mesh->med[mesh->type[t]].mua;
	musp=mesh->med[mesh->type[t]].mus;
	nref=mesh->med[mesh->type[t]].n;

	dcoeff=1.0/(3.0*(mua+musp));

        /*calculate del_x,del_y and del_z for 4 basis functions per elem*/
        derx[0]=-((N[i3].y*N[i4].z-N[i3].z*N[i4].y)-N[i2].y*(N[i4].z-N[i3].z)
               +N[i2].z*(N[i4].y-N[i3].y))*RVe6;
        dery[0]=((N[i3].x*N[i4].z-N[i4].x*N[i3].z)-N[i2].x*(N[i4].z-N[i3].z)
               +N[i2].z*(N[i4].x-N[i3].x))*RVe6;
        derz[0]=-((N[i3].x*N[i4].y-N[i3].y*N[i4].x)-N[i2].x*(N[i4].y-N[i3].y)
               +N[i2].y*(N[i4].x-N[i3].x))*RVe6;

        derx[1]=((N[i3].y*N[i4].z-N[i3].z*N[i4].y)-N[i1].y*(N[i4].z-N[i3].z)
               +N[i1].z*(N[i4].y-N[i3].y))*RVe6;
        dery[1]=-((N[i3].x*N[i4].z-N[i4].x*N[i3].z)-N[i1].x*(N[i4].z-N[i3].z)
               +N[i1].z*(N[i4].x-N[i3].x))*RVe6;
        derz[1]=((N[i3].x*N[i4].y-N[i3].y*N[i4].x)-N[i1].x*(N[i4].y-N[i3].y)
               +N[i1].y*(N[i4].x-N[i3].x))*RVe6;

        derx[2]=-((N[i2].y*N[i4].z-N[i2].z*N[i4].y)-N[i1].y*(N[i4].z-N[i2].z)
               +N[i1].z*(N[i4].y-N[i2].y))*RVe6;
        dery[2]=((N[i2].x*N[i4].z-N[i4].x*N[i2].z)-N[i1].x*(N[i4].z-N[i2].z)
               +N[i1].z*(N[i4].x-N[i2].x))*RVe6;
        derz[2]=-((N[i2].x*N[i4].y-N[i2].y*N[i4].x)-N[i1].x*(N[i4].y-N[i2].y)
               +N[i1].y*(N[i4].x-N[i2].x))*RVe6;

        derx[3]=((N[i2].y*N[i3].z-N[i2].z*N[i3].y)-N[i1].y*(N[i3].z-N[i2].z)
               +N[i1].z*(N[i3].y-N[i2].y))*RVe6;
        dery[3]=-((N[i2].x*N[i3].z-N[i3].x*N[i2].z)-N[i1].x*(N[i3].z-N[i2].z)
               +N[i1].z*(N[i3].x-N[i2].x))*RVe6;
        derz[3]=((N[i2].x*N[i3].y-N[i2].y*N[i3].x)-N[i1].x*(N[i3].y-N[i2].y)
               +N[i1].y*(N[i3].x-N[i2].x))*RVe6;

        /*loop over index i*/
        for(i=0;i<4;i++){
          ii = ((int *)(mesh->elem+t))[i];
          /*loop over index j*/
          for(j=i;j<4;j++){
            jj =((int *)(mesh->elem+t))[j];
            deldotdel=(derx[i]*derx[j]+dery[i]*dery[j]+derz[i]*derz[j])*Ve;
	    
	    if(fem->deldotdel)
	        fem->deldotdel[t*10+(i<<2)+j-i-(i>1)-((i>2)<<1)]=deldotdel;

            sm = Ve*0.05;
            if (i==j) sm = Ve*0.1;

            /* D*Ve*delFi*delF */
            ra=dcoeff*deldotdel+mua*sm;
	    ri=cfg->omega*R_C0*nref*sm;

            /*add the values to the vectorized sparse matrix*/
            if(ii==jj){
                fem->Dr[ii]+=ra;
		if(fem->Di)
		    fem->Di[ii]+=ri;
            }else{
                ij=sub2ind(mesh,ii,jj);
                fem->Ar[ij]+=ra;
		if(fem->Ai)
		    fem->Ai[ij]+=ri;		    
            }
          }
        }
    }
}

int sub2ind(tetmesh *mesh, int i, int j){
    int base,ij=-1,k;
    base = (int)(mesh->idxsum[i] - mesh->idxcount[i]);
    for(k=0;k<(int)(mesh->idxcount[i]);k++){
         ij = k + base;
         if ((int)(mesh->cols[ij]) == j) break;
    }
    return ij;
}

extern "C" int rb3_throw_exception(const int id, const char *msg, const char *filename, const int linenum){
     printf("Redbird ERROR (%d): %s in unit %s:%d\n",id,msg,filename,linenum);
     throw(msg);
     return id;
}

void rb3_usage(){
     printf("Usage:\n    [flux,detphoton]=rbforward(cfg);\n\nPlease run 'help mmclab' for more details.\n");
}

extern "C" void rb3_flush(){
#ifndef MATLAB_MEX_FILE
	mexEvalString("fflush(stdout);");
#else
	mexEvalString("pause(.0001);");
#endif
}
