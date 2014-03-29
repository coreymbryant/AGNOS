


namespace AGNOS{
  /** Comparison functor for index sets */
  bool indexSetCompare( 
      const std::vector<unsigned int>& a, 
      const std::vector<unsigned int>& b 
      )
  {
    assert( a.size() == b.size() ) ;
    bool aLessThanb = false;

    for(unsigned int i=0; i<a.size(); i++)
      if ( a[i] < b[i] )
      {
        aLessThanb = true ;
        break;
      }
      else if ( a[i] > b[i] )
        break;

    return aLessThanb;
  }

}


#endif // AGNOS_DEFINES_H
