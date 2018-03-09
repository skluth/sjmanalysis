#ifndef ANALYSIS_HH
#define ANALYSIS_HH

#include <string>

class Analysis {
public:
  Analysis( const std::string & src,
	    const std::string & rec,
	    const std::string & cts,
	    const std::string & mcc="none",
	    const std::string & r2="none",
	    const std::string & bkgsts="none",
	    const std::string & unfsrc="none",
	    const std::string & unfm="none" );
  Analysis( const std::string & );
  Analysis();
  ~Analysis();
  std::string getSource() const;
  std::string getReco() const;
  std::string getCuts() const;
  std::string getMccuts() const;
  std::string getReco2() const;
  std::string getBkgStatus() const;
  std::string getUnfoldSource() const;
  std::string getUnfoldMethod() const;
  std::string getTag() const;
  void Print() const;
  void setSource( const std::string & );
  void setReco( const std::string & );
  void setCuts( const std::string & );
  void setMccuts( const std::string & );
  void setReco2( const std::string & );
  void setBkgStatus( const std::string & );
  void setUnfoldSource( const std::string & );
  void setUnfoldMethod( const std::string & );
private:
  std::string source;
  std::string reco;
  std::string cuts;
  std::string mccuts;
  std::string reco2;
  std::string bkgstatus;
  std::string unfoldsource;
  std::string unfoldmethod;
};

#endif
