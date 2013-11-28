#include "mex.h"
#include "clipper.hpp"

using namespace ClipperLib;
using namespace std;

double PsArea(Polygons & ps){
  double resarea=0;
  for(int i = 0; i<ps.size(); i++){
    resarea+=Area(ps[i]);
  }
  return resarea<0?-resarea:resarea;
}

Polygons* genPoly(long64* rects, int starti, int fini,int len){
  Polygons* solution = NULL;
  Polygons* solution_old = NULL;
  //mexPrintf("1\n");
  double prevtotal=0,total=0;
  for(int i = starti; i < fini; i++){
    //mexPrintf("2\n");
    Polygons* rect = new Polygons();
    rect->resize(1);
    (*rect)[0].resize(4);
    //mexPrintf("4\n");
    (*rect)[0][0].X=(long64)rects[i];
    (*rect)[0][0].Y=(long64)rects[i+1*len];
    (*rect)[0][1].X=(long64)rects[i];
    (*rect)[0][1].Y=(long64)rects[i+3*len];
    (*rect)[0][2].X=(long64)rects[i+2*len];
    (*rect)[0][2].Y=(long64)rects[i+3*len];
    (*rect)[0][3].X=(long64)rects[i+2*len];
    (*rect)[0][3].Y=(long64)rects[i+1*len];
    //mexPrintf("3\n");
    if(i==starti){
      solution_old=rect;
    }else{
      Clipper c;
      c.AddPolygons(*solution_old,ptSubject);
      c.AddPolygons(*rect,ptClip);
      solution = new Polygons();
      c.Execute(ctUnion,*solution,pftNonZero,pftNonZero);
      //mexPrintf("arean %f\n",total);
      delete solution_old;
      delete rect;
      solution_old=solution;
    }
  }
  return solution_old;
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
  mexPrintf("entered mexfn\n");
  mexEvalString("drawnow;");
  long64* rects=((long64*)mxGetPr(prhs[0]));
  int len = *(mxGetDimensions(prhs[0]));
  long64* ntofind=((long64*)mxGetPr(prhs[1]));
  //mexPrintf("len %d\n",len);
  int out[1];
  out[0] = (int) (*ntofind);
          mexPrintf("outsize:%d\n",out[0]);
          mexEvalString("drawnow;");
  plhs[0]=mxCreateNumericArray(1, out, mxDOUBLE_CLASS, mxREAL);
  plhs[1]=mxCreateNumericArray(1, out, mxDOUBLE_CLASS, mxREAL);
  double* outdata=((double*)mxGetPr(plhs[0]));
  double* outscores=((double*)mxGetPr(plhs[1]));
  int maximg=0;
  int maxdetr=0;
  for(int i=0; i<len; i++){
    if(rects[i+4*len]>maxdetr) maxdetr=(long64)rects[i+4*len];
    if(rects[i+5*len]>maximg) maximg=(long64)rects[i+5*len];
  }
  vector<vector<Polygons*>*> allpolys(maxdetr+1);
  for(int i = 0; i<=maxdetr; i++){
    allpolys[i] = new vector<Polygons*>(maximg+1);
  }
  int i = 0;
  mexPrintf("generating polys_new\n");
  mexEvalString("drawnow;");
  while(i<len){
    int currdetr=(long64)rects[i+4*len];
    int currimg=(long64)rects[i+5*len];
    int j = i+1;
    while(j<len && (long64)rects[j+4*len]==currdetr && (long64)rects[j+5*len]==currimg){
      j++;
    }
    (*allpolys[currdetr])[currimg]=genPoly(rects,i,j,len);
    i=j;
  }
  vector<Polygons*> totalpolys(maximg+1);
  vector<bool> tprequiresdelete(maximg+1,false);
  vector<double> prevcontribution(maxdetr+1,1.0/0.0);
  mexPrintf("computing unions\n");
  mexEvalString("drawnow;");
  double totalcoverage=0;
  for(int i = 0; i<*ntofind; ++i){
    int bestdetr=0;
    double bestarea=0;
    vector<Polygons*> bestunionpolys(maximg+1);
    vector<bool> buprequiresdelete(maximg+1,false);
    vector<bool> bupcamefromtp(maximg+1,false);
    for(int curdetr=1; curdetr<=maxdetr; curdetr++){
      if(prevcontribution[curdetr]<bestarea-totalcoverage){
        continue;
      }
          //mexPrintf("curdetr:%d\n",curdetr);
          //mexEvalString("drawnow;");
      double totalarea=0;
      vector<Polygons*> tmpunionpolys(maximg+1);
      vector<bool> tuprequiresdelete(maximg+1,false);
      vector<bool> tupcamefromtp(maximg+1,false);
      for(int curimg=1; curimg<=maximg; curimg++){
        Clipper c;
        vector<Polygons*>* curdetrpolys;
        curdetrpolys=(allpolys[curdetr]);
        if(totalpolys[curimg]==NULL && (*curdetrpolys)[curimg]==NULL){
        }else if(totalpolys[curimg]==NULL && (*curdetrpolys)[curimg]!=NULL){
          //if(i>0){
          //mexPrintf("from new\n");
          //mexEvalString("drawnow;");
          //}
          tmpunionpolys[curimg]=(*curdetrpolys)[curimg];
          totalarea+=PsArea(*(*curdetrpolys)[curimg]);
        }else if(totalpolys[curimg]!=NULL && (*curdetrpolys)[curimg]==NULL){
          //mexPrintf("from prev\n");
          //mexEvalString("drawnow;");
          tmpunionpolys[curimg]=totalpolys[curimg];
          totalarea+=PsArea(*(totalpolys[curimg]));
          tuprequiresdelete[curimg]=tprequiresdelete[curimg];
          tupcamefromtp[curimg]=true;
        }else{
          //mexPrintf("unioning\n");
          //mexEvalString("drawnow;");
          c.AddPolygons(*(totalpolys[curimg]),ptSubject);
          c.AddPolygons(*((*curdetrpolys)[curimg]),ptClip);
          tmpunionpolys[curimg]=new Polygons();
          c.Execute(ctUnion,*(tmpunionpolys[curimg]),pftNonZero,pftNonZero);
          tuprequiresdelete[curimg]=true;
          totalarea+=PsArea(*tmpunionpolys[curimg]);
          //mexPrintf("doneunion\n");
          //mexEvalString("drawnow;");
        }
      }
      prevcontribution[curdetr]=totalarea-totalcoverage;
      if(totalarea>bestarea){
        bestarea=totalarea;
        for(int curimg=1; curimg<=maximg; curimg++){
          if(buprequiresdelete[curimg] && !bupcamefromtp[curimg]) delete bestunionpolys[curimg];
        }
        bestunionpolys=tmpunionpolys;
        buprequiresdelete=tuprequiresdelete;
        bupcamefromtp=tupcamefromtp;
        bestdetr=curdetr;
        //mexPrintf("new best: %d, area %f\n", bestdetr,bestarea);
        //mexEvalString("drawnow;");
      }else{
        for(int curimg=1; curimg<=maximg; curimg++){
          if(tuprequiresdelete[curimg] && !tupcamefromtp[curimg]) delete tmpunionpolys[curimg];
        }
      }
    }
    outdata[i]=bestdetr;
    mexPrintf("selected: %d, area %f; completed %d / %d\n", bestdetr, bestarea, i+1,*ntofind);
    mexEvalString("drawnow;");
    outscores[i]=bestarea-totalcoverage;
    totalcoverage=bestarea;
    for(int curimg=1; curimg<=maximg; curimg++){
      if(tprequiresdelete[curimg] && !bupcamefromtp[curimg]) delete totalpolys[curimg];
    }
    totalpolys=bestunionpolys;
    tprequiresdelete=buprequiresdelete;
  }
  for(int curimg=1; curimg<=maximg; curimg++){
    if(tprequiresdelete[curimg]) delete totalpolys[curimg];
    for(int curdetr=1; curdetr<=maxdetr; curdetr++){
      //mexPrintf("allpolydelete:\n");
      //mexEvalString("drawnow;");
      if((*allpolys[curdetr])[curimg]!=NULL) delete (*allpolys[curdetr])[curimg];
    }
    //  mexPrintf("allpolyfinaldelete:\n");
    //  mexEvalString("drawnow;");
  }
  for(int curdetr=1;curdetr<=maxdetr;curdetr++){
    delete allpolys[curdetr];
  }

}

