#ifndef ANALYSIS_HH
#define ANALYSIS_HH

#include <string>
using std::string;

class Analysis {
public:
  Analysis( const string& source, const string& reco, const string& cuts, 
	    const string& mccuts="none", const string& reco2="none",
	    const string& unfoldsource="none", const string& unfoldmethod="none" );
  Analysis();
  ~Analysis();
  string getSource() const;
  string getReco() const;
  string getCuts() const;
  string getMccuts() const;
  string getReco2() const;
  string getUnfoldSource() const;
  string getUnfoldMethod() const;
  string getTag() const;
  void print() const;
private:
  string a_source;
  string a_reco;
  string a_cuts;
  string a_mccuts;
  string a_reco2;
  string a_unfoldsource;
  string a_unfoldmethod;
};

#endif
