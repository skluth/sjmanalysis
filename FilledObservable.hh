#ifndef FILLEDOBSERVABLE_HH
#define FILLEDOBSERVABLE_HH

#include "Rtypes.h"
#include "Analysis.hh"

#include <string>
#include <map>

class DataStructure;
class DifferentialDataStructure;
class JetrateDataStructure;
class MatrixDataStructure;


class FilledObservable {

public:

  FilledObservable( const std::string &, 
  		    const std::map<std::string,DifferentialDataStructure*> &,
  		    const std::map<std::string,MatrixDataStructure*> & mds= 
  		    std::map<std::string,MatrixDataStructure*>() );

  FilledObservable( const std::string &, 
  		    const std::map<std::string,JetrateDataStructure*> & );

  ~FilledObservable() {}
  void finalise();
  void Print() const;
  const std::map<std::string,DataStructure*> & getData() const;
  const std::map<std::string,MatrixDataStructure*> & getMigrationMatrices() const;
  std::map<std::string,MatrixDataStructure*> getErrorMatrices() const;
  DataStructure* getDataStructure( const Analysis & ) const;
  void setDataStructure( DataStructure*, const Analysis &  );
  std::string getName() const { return name; }
  bool containsAnalysis( const Analysis & ) const;
  void requestErrorMatrices() { lerrorMatrices= true; }
  bool makeErrorMatrices() const { return lerrorMatrices; }
  
private:

  std::string name;
  bool lerrorMatrices;
  std::map<std::string,DataStructure*> datastructures;
  std::map<std::string,MatrixDataStructure*> migrationMatrices;

};

#endif
