#include <bits/stdc++.h>
using namespace std;
using Chain=vector<uint16_t>;
thread_local int N,lim,a[24]; thread_local vector<Chain> out;
void save(){Chain v; for(int i=0;i<lim;i++)v.push_back(a[i]);v.push_back(N);out.push_back(move(v));}
bool sumexists(int j,int x){int l=0,r=j;while(l<=r){int s=a[l]+a[r];if(s==x)return true;if(s<x)l++;else r--;}return false;}
void dfs(int j){int r=lim-j, m=a[j]; if((m<<r)<N)return;
 if(r==1){if(sumexists(j,N))save();return;}
 if(r==2){ // for shortest chains, the last new term must participate in N
   int minc=(N+1)/2;
   for(int k=0;k<=j;k++){int c=N-a[k];if(c<=m||c>2*m)continue;if(sumexists(j,c)){a[j+1]=c;save();}}
   if(N%2==0 && N/2>m && N/2<=2*m && sumexists(j,N/2)){a[j+1]=N/2;save();}
   return;
 }
 int low=max(m+1,(N+(1<<(r-1))-1)>>(r-1)),high=N-1;
 // local unique candidates via bitset
 uint64_t seen[24]={};
 for(int i=j;i>=0;i--){if(a[i]*2<low)break;for(int k=i;k>=0;k--){int x=a[i]+a[k];if(x<low)break;if(x>high)continue;auto &w=seen[x>>6];uint64_t b=1ULL<<(x&63);if(w&b)continue;w|=b;a[j+1]=x;dfs(j+1);}}
}
vector<Chain> generate(int n){N=n;a[0]=1;for(lim=0;(1<<lim)<N;lim++); for(;;lim++){out.clear();dfs(0);if(!out.empty())return move(out);}}
bool compatible(const Chain &a,const Chain &b){ // shared prefix values; for every shared length require a pair in intersection
 int common[24], nc=0,i=0,j=0;
 while(i<(int)a.size()&&j<(int)b.size()){
 if(a[i]<b[j]){i++;continue;} if(a[i]>b[j]){j++;continue;}
 int x=a[i];if(x!=1){int l=0,r=nc-1;bool ok=false;while(l<=r){int s=common[l]+common[r];if(s==x){ok=true;break;}if(s<x)l++;else r--;}if(!ok)return false;}
 common[nc++]=x;i++;j++;
 }return true;
}
int main(int argc,char**argv){if(argc>1){for(int i=1;i<argc;i++){int n=atoi(argv[i]);auto t=chrono::steady_clock::now();auto v=generate(n);cerr<<n<<" ASI "<<v[0].size()-1<<" chains "<<v.size()<<" sec "<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<endl; if(v.size()<20)for(auto &c:v){for(auto x:c)cerr<<x<<",";cerr<<endl;}}return 0;}
 int K=1387;vector<vector<Chain>> plans(K+1); auto start=chrono::steady_clock::now();
 #pragma omp parallel for schedule(dynamic,1)
 for(int n=2;n<=K;n++){plans[n]=generate(n);if(n%100==0){
 #pragma omp critical
 cerr<<"generated "<<n<<" ("<<plans[n].size()<<") elapsed "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;}}
 size_t total=0;for(int n=2;n<=K;n++)total+=plans[n].size();cerr<<"generated all "<<total<<" chains in "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<" seconds\n";
 long long tested=0,checks=0; for(int n=3;n<=K;n++){for(int m=2;m<n;m++){bool ok=false;auto &x=plans[m],&y=plans[n]; for(auto &a:x){for(auto &b:y){checks++;if(compatible(a,b)){ok=true;break;}}if(ok)break;}tested++;if(!ok)cout<<"NO SOJAS "<<m<<" "<<n<<" chains "<<x.size()<<" "<<y.size()<<" pairs tested "<<tested<<"\n"<<flush;}if(n%100==0||n==1386||n==1387)cerr<<"checked through "<<n<<" pairs "<<tested<<" compatibility tests "<<checks<<" elapsed "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;}
}