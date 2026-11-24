/* fclq.c - Fable review 2026-11-24. Independent clique / independent-set decision.
   Method: branch and bound with greedy colour classes (Tomita style): at a node needing r more
   vertices from candidate set P, r-1 greedy colour classes (independent sets) are peeled off P;
   any r-clique in P must use a vertex outside them, so only the remaining vertices are branched on.
   Usage: check_branch_bound_fable FILE T c|i [threads]    FILE = n rows of n chars 0/1.  c: clique of size T?  i: independent set of size T?
   Prints nodes, and a verified witness if one exists. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#define MAXW 32
typedef uint64_t u64;
static int n, W;
static u64 (*adj)[MAXW];
static volatile int found = 0;
static int witness[64];
static inline int popc(const u64 *a){int c=0;for(int i=0;i<W;i++)c+=__builtin_popcountll(a[i]);return c;}
static inline int first(const u64 *a){for(int i=0;i<W;i++)if(a[i])return i*64+__builtin_ctzll(a[i]);return -1;}
/* peel r-1 greedy colour classes from P; write the rest to B */
static void peel(const u64 *P, int r, u64 *B){
  u64 U[MAXW], Q[MAXW]; memcpy(U,P,W*8);
  for(int c=0;c<r-1;c++){
    memcpy(Q,U,W*8); int v;
    int any=0;
    while((v=first(Q))>=0){ any=1; U[v>>6]&=~(1ULL<<(v&63)); Q[v>>6]&=~(1ULL<<(v&63)); for(int i=0;i<W;i++)Q[i]&=~adj[v][i]; }
    if(!any)break;
  }
  memcpy(B,U,W*8);
}
static int expand(u64 *P, int r, int *stack, int depth, u64 *nodes){
  (*nodes)++;
  if(found) return 0;
  if(r==0) return 1;
  int cnt=popc(P); if(cnt<r) return 0;
  if(r==1){ stack[depth]=first(P); return 1; }
  if(r==2){
    u64 Q[MAXW]; memcpy(Q,P,W*8); int v;
    while((v=first(Q))>=0){ Q[v>>6]&=~(1ULL<<(v&63));
      for(int i=0;i<W;i++){ u64 x=Q[i]&adj[v][i]; if(x){ stack[depth]=v; stack[depth+1]=i*64+__builtin_ctzll(x); return 1; } } }
    return 0;
  }
  u64 B[MAXW], N[MAXW]; peel(P,r,B);
  int v;
  while((v=first(B))>=0){
    B[v>>6]&=~(1ULL<<(v&63));
    for(int i=0;i<W;i++)N[i]=P[i]&adj[v][i];
    stack[depth]=v;
    if(expand(N,r-1,stack,depth+1,nodes)) return 1;
    P[v>>6]&=~(1ULL<<(v&63));
    if(found) return 0;
  }
  return 0;
}
int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage: check_branch_bound_fable FILE T c|i [threads]\n");return 2;}
  int T=atoi(argv[2]); char mode=argv[3][0]; int th=argc>4?atoi(argv[4]):4;
  FILE*f=fopen(argv[1],"r"); if(!f){perror("open");return 2;}
  char *line=malloc(1<<16); int rows=0; char **M=NULL;
  while(fgets(line,1<<16,f)){ int L=strlen(line); while(L>0 ? (line[L-1]=='\n'||line[L-1]=='\r') : 0)line[--L]=0; {int k=0;for(int q=0;q<L;q++)if(line[q]!=' ')line[k++]=line[q];line[k]=0;L=k;} if(!L)continue;
    M=realloc(M,(rows+1)*sizeof(char*)); M[rows]=strdup(line); rows++; }
  fclose(f); n=rows; W=(n+63)/64; if(W>MAXW){fprintf(stderr,"too big\n");return 2;}
  long edges=0;
  for(int i=0;i<n;i++){ if((int)strlen(M[i])!=n){fprintf(stderr,"row %d length %zu != %d\n",i,strlen(M[i]),n);return 2;}
    if(M[i][i]!='0'){fprintf(stderr,"nonzero diagonal at %d\n",i);return 2;}
    for(int j=0;j<n;j++){ if(M[i][j]!='0'){ if(M[i][j]!='1'){fprintf(stderr,"bad char\n");return 2;} } if(M[i][j]!=M[j][i]){fprintf(stderr,"not symmetric %d %d\n",i,j);return 2;} if(M[i][j]=='1')edges++; } }
  edges/=2;
  adj=calloc(n,sizeof(*adj));
  for(int i=0;i<n;i++)for(int j=0;j<n;j++) if(i!=j){ int e=(M[i][j]=='1'); if(mode=='i')e=!e; if(e)adj[i][j>>6]|=1ULL<<(j&63); }
  printf("file=%s n=%d edges=%ld symmetric=yes zero_diag=yes mode=%s T=%d threads=%d\n",argv[1],n,edges,mode=='i'?"independent":"clique",T,th);
  fflush(stdout);
  u64 P[MAXW]; memset(P,0,sizeof P); for(int i=0;i<n;i++)P[i>>6]|=1ULL<<(i&63);
  u64 B[MAXW]; peel(P,T,B);
  int nb=popc(B); int *bv=malloc(nb*sizeof(int)); u64 (*PP)[MAXW]=calloc(nb,sizeof(*PP));
  { int v,k=0; while((v=first(B))>=0){ B[v>>6]&=~(1ULL<<(v&63)); bv[k]=v; for(int i=0;i<W;i++)PP[k][i]=P[i]&adj[v][i]; P[v>>6]&=~(1ULL<<(v&63)); k++; } }
  printf("root: %d branch vertices after peeling %d colour classes\n",nb,T-1); fflush(stdout);
  u64 total=0; double t0=omp_get_wtime();
  omp_set_num_threads(th);
  #pragma omp parallel for schedule(dynamic,1) reduction(+:total)
  for(int k=0;k<nb;k++){
    if(found) continue;
    int stack[64]; u64 nodes=0; u64 N[MAXW]; memcpy(N,PP[k],W*8);
    stack[0]=bv[k];
    if(expand(N,T-1,stack,1,&nodes)){
      #pragma omp critical
      { if(!found){ found=1; memcpy(witness,stack,T*sizeof(int)); } }
    }
    total+=nodes;
  }
  printf("nodes=%llu wall=%.1fs\n",(unsigned long long)total,omp_get_wtime()-t0);
  if(found){
    int ok=1; for(int a=0;a<T;a++)for(int b=a+1;b<T;b++){ int x=witness[a],y=witness[b]; if(x==y||!((adj[x][y>>6]>>(y&63))&1))ok=0; }
    printf("witness:"); for(int a=0;a<T;a++)printf(" %d",witness[a]); printf("  verified=%s\n",ok?"yes":"NO");
    printf("RESULT: %s of size %d exists: YES\n",mode=='i'?"independent set":"clique",T);
  } else printf("RESULT: %s of size %d exists: no\n",mode=='i'?"independent set":"clique",T);
  return 0;
}
