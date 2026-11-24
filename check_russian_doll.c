/* Independent checker written by the main session (2026-11-23), not derived from check.py/check2.py/check3.py/indep.c.
   Reads a 0/1 adjacency matrix (any non-digit chars ignored), checks symmetry and zero diagonal,
   counts K4 by edge + common-neighbour bitsets, and decides whether an independent set of size T exists
   using Ostergard's Russian-doll algorithm (c[i] = largest independent set inside vertices i..n-1). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#define MAXN 1600
#define W ((MAXN+63)/64)
typedef struct { uint64_t w[W]; } bs;
static int n, nw, T;
static bs adj[MAXN], nadj[MAXN];
static int c[MAXN];
static int found;
static long long nodes;
static int pc(const bs *a){int s=0;for(int i=0;i<nw;i++)s+=__builtin_popcountll(a->w[i]);return s;}
static int lowest(const bs *a){for(int i=0;i<nw;i++)if(a->w[i])return i*64+__builtin_ctzll(a->w[i]);return -1;}
/* cand = vertices still allowed (all > current doll start), size = chosen so far */
static void expand(bs *cand, int size, int target){
  nodes++;
  int k = pc(cand);
  if(k==0){ if(size>=target) found=1; return; }
  while(k>0){
    if(size + k < target) return;
    int v = lowest(cand);
    if(size + c[v] < target) return;
    cand->w[v/64] &= ~(1ULL<<(v%64));
    bs nc; for(int i=0;i<nw;i++) nc.w[i] = cand->w[i] & nadj[v].w[i];
    expand(&nc, size+1, target);
    if(found) return;
    k--;
  }
}
int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage: check_russian_doll graph.txt T c|i\n");return 2;}
  T=atoi(argv[2]);
  FILE*f=fopen(argv[1],"r"); if(!f){perror("open");return 2;}
  static int vals[MAXN*MAXN]; long m=0; int ch;
  while((ch=fgetc(f))!=EOF){ if(ch=='0'||ch=='1'){ if(m>=(long)MAXN*MAXN){fprintf(stderr,"too big\n");return 2;} vals[m++]=ch-'0'; } }
  fclose(f);
  n=0; while((long)n*n<m) n++;
  if((long)n*n!=m){fprintf(stderr,"not square: %ld entries\n",m);return 2;}
  nw=(n+63)/64;
  int bad=0; long edges=0;
  for(int i=0;i<n;i++)for(int j=0;j<n;j++){
    int a=vals[(long)i*n+j];
    if(i==j&&a){bad=1;} if(a!=vals[(long)j*n+i]) bad=1;
    if(a&&i!=j){adj[i].w[j/64]|=1ULL<<(j%64); if(i<j)edges++;}
    if(!a&&i!=j){nadj[i].w[j/64]|=1ULL<<(j%64);}
  }
  printf("n=%d edges=%ld symmetric_and_zero_diag=%s\n",n,edges,bad?"NO":"yes");
  if(bad) return 1;
  int cmode = (argv[3][0]=='c');
  const char *what = cmode ? "clique number" : "independence number";
  if(cmode){ for(int i=0;i<n;i++){bs t=adj[i];adj[i]=nadj[i];nadj[i]=t;} printf("MODE clique: searching for a clique of size T (complement independence)\n"); } else printf("MODE independent set of size T\n");
  /* K4: for each edge u<v, common neighbours w>v, need an edge inside them */
  long k4=0;
  if(n<=300)
  for(int u=0;u<n;u++)for(int v=u+1;v<n;v++) if(adj[u].w[v/64]>>(v%64)&1){
    bs cm; for(int i=0;i<nw;i++) cm.w[i]=adj[u].w[i]&adj[v].w[i];
    for(int w=v+1;w<n;w++) if(cm.w[w/64]>>(w%64)&1){
      bs x; for(int i=0;i<nw;i++) x.w[i]=cm.w[i]&adj[w].w[i];
      for(int z=w+1;z<n;z++) if(x.w[z/64]>>(z%64)&1) k4++;
    }
  }
  if(n<=300) printf("K4 count: %ld\n",k4); else printf("K4 count: skipped (n > 300)\n");
  clock_t t0=clock();
  /* Russian doll: c[i] is computed for i = n-1 down to 0, with the exact max in vertices i..n-1 */
  int best=0;
  for(int i=n-1;i>=0;i--){
    bs cand; memset(&cand,0,sizeof cand);
    for(int j=i+1;j<n;j++) if(nadj[i].w[j/64]>>(j%64)&1) cand.w[j/64]|=1ULL<<(j%64);
    /* try to improve best by 1 using vertex i */
    found=0;
    expand(&cand,1,best+1);
    if(found) best++;
    c[i]=best;
    if(best>=T) break;
  }
  printf("%s %s %d (target %d)  nodes=%lld  cpu=%.1fs\n", what, best>=T?">=":"=", best, T, nodes, (double)(clock()-t0)/CLOCKS_PER_SEC);
  printf("%s of size %d exists: %s\n", cmode?"clique":"independent set", T, best>=T?"YES":"no");
  return 0;
}
