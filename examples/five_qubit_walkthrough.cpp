// Walk-through for docs/dp-image.tex, section 2: the [[5,1,3]] code and a 5-element error set.
// Prints the diagram of the set, the image I(u) of every node (the DP by hand, on the real
// diagram) and the verdicts of both methods.
#include <spbdd/spbdd.hpp>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <vector>
using namespace spbdd;
typedef std::set<unsigned> VS;           // set of 6-bit vectors: bits 0-3 sigma, bit 4,5 lambda
static unsigned W[10];
static std::map<int,int> id; static std::vector<Bdd> nodes;
static int nid(const Bdd&u){ if(u.is_false())return -1; if(u.is_true())return -2; auto it=id.find(u.node()); if(it!=id.end())return it->second; int k=nodes.size(); id[u.node()]=k; nodes.push_back(u); nid(u.low()); nid(u.high()); return k; }
static std::string vs(unsigned v){ std::string s; for(int i=0;i<4;i++) s+=(v>>i&1)?'1':'0'; s+='|'; for(int i=4;i<6;i++) s+=(v>>i&1)?'1':'0'; return s; }
static std::string show(const VS&s){ std::string o="{"; bool f=true; for(unsigned v:s){ if(!f)o+=", "; f=false; o+=vs(v);} return o+"}"; }
static std::map<std::pair<int,int>,VS> memo;
static Manager* M;
static VS img(const Bdd&u,int lvl){
  if(u.is_false()) return {};
  if(lvl==10) return {0u};
  int key=u.is_true()?-2:id[u.node()];
  auto it=memo.find({key,lvl}); if(it!=memo.end()) return it->second;
  int v=M->level_to_var(lvl); VS out;
  if(!u.is_constant() && M->var_to_level(u.top_var())==lvl){ VS lo=img(u.low(),lvl+1), hi=img(u.high(),lvl+1); out=lo; for(unsigned x:hi) out.insert(x^W[v]); }
  else { VS h=img(u,lvl+1); out=h; for(unsigned x:h) out.insert(x^W[v]); }
  memo[{key,lvl}]=out; return out; }
int main(int argc,char**argv){
  PauliSpace sp(5); M=&sp.manager();
  StabilizerCode code(sp,{"XZZXI","IXZZX","XIXZZ","ZXIXZ"});
  for(int q=0;q<5;q++)for(int b=0;b<2;b++){std::string s(5,'I');s[q]=b?'Z':'X';auto sy=code.syndrome(s);auto sg=code.logical_signature(s);unsigned p=0;for(int i=0;i<4;i++)if(sy[i])p|=1u<<i;for(int j=0;j<2;j++)if(sg[j])p|=1u<<(4+j);W[2*q+b]=p;}
  std::vector<std::string> L={"IIIII","IIXII","IIIYX","XIIII","IZIII"};
  PauliSet E=sp.from(L);
  printf("|E|=%.0f nodes=%zu\n",E.size(),E.node_count());
  Bdd root=E.bdd(); int r=nid(root);
  // print nodes
  for(size_t k=0;k<nodes.size();k++){ const Bdd&u=nodes[k]; int v=u.top_var(); auto nm=[&](int i){return i==-1?std::string("F"):i==-2?std::string("T"):"n"+std::to_string(i);};
    printf("n%zu: var %d (%c%d, level %d)  low=%s high=%s\n",k,v,v%2?'Z':'X',v/2+1,M->var_to_level(v),nm(nid(u.low())).c_str(),nm(nid(u.high())).c_str()); }
  printf("root=n%d\n",r);
  // image per node at its own level
  for(size_t k=0;k<nodes.size();k++){ int lv=M->var_to_level(nodes[k].top_var()); printf("I(n%zu,level %d) = %s\n",k,lv,show(img(nodes[k],lv)).c_str()); }
  VS G=img(root,0); printf("G = %s\n",show(G).c_str());
  auto p=code.find_inequivalent_pair(E); printf("library DP: %s", p?"pair":"no pair"); if(p) printf("  %s , %s\n",p->first.c_str(),p->second.c_str()); else printf("\n");
  auto q=code.find_inequivalent_pair(E,StabilizerCode::PairMethod::Compose); printf("library compose: %s\n", q?"pair":"no pair");
  for(auto&s:L){auto sy=code.syndrome(s);auto sg=code.logical_signature(s);std::string a,b;for(bool x:sy)a+=x?'1':'0';for(bool x:sg)b+=x?'1':'0';printf("%s sigma=%s lambda=%s\n",s.c_str(),a.c_str(),b.c_str());}
}
