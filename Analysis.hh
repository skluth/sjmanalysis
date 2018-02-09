#ifndef ANALYSIS_HH
#define ANALYSIS_HH

#include <string>
using std::string;

class Analysis {
public:
  Analysis( const string&, const string&, const string&, 
	    const string& mcc="none", const string& r2="none",
	    const string& bkgsts="none",
	    const string& unfsrc="none", const string& unfm="none" );
  Analysis( const string& );
  Analysis();
  ~Analysis();
  string getSource() const;
  string getReco() const;
  string getCuts() const;
  string getMccuts() const;
  string getReco2() const;
  string getBkgStatus() const;
  string getUnfoldSource() const;
  string getUnfoldMethod() const;
  string getTag() const;
  void print() const;
  void setSource( const string& );
  void setReco( const string& );
  void setCuts( const string& );
  void setMccuts( const string& );
  void setReco2( const string& );
  void setBkgStatus( const string& );
  void setUnfoldSource( const string& );
  void setUnfoldMethod( const string& );
private:
  string source;
  string reco;
  string cuts;
  string mccuts;
  string reco2;
  string bkgstatus;
  string unfoldsource;
  string unfoldmethod;
};

#endif
