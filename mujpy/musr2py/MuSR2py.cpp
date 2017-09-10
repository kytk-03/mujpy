// Test musr2py
#include <new> //For std::nothrow
#include <inttypes.h> //For copying memory
#include "MuSR_td_PSI_bin.h"

extern "C"  //Tells the compile to use C-linkage for the next scope.
{
    void * init( void )
    {
        // New  MuSR_td_PSI_bin class
        return new(std::nothrow) MuSR_td_PSI_bin; // puntatore a finta classe c che si wrappa sulla classe c++
    }

    void destroy (void *ptr)
    {
        // Avoid throwing exceptions. How?
        delete (MuSR_td_PSI_bin*) ptr;
    }
    
    
    // Read File
    int read(void *ptr, char * filename)  // void*ptr è l'inidizzo di memoria dove risiede la classe di Amato 
                                          //  (utilizzato sempre con MuSR_td_PSI_bin della nuova classe c una decina di righe qui sopra 
    {
        // Avoid throwing exceptions.  
        try 
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            return ref->read(filename);
        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }
    
    
    // get length of histogram in bins
    int get_histoLength_bin(void *ptr)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            return ref->get_histoLength_bin();
        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }

    // get number of histograms 
    int get_numberHisto_int(void *ptr)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            return ref->get_numberHisto_int();
        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }

    // get histogram as array of integers    
    int get_histo_array_int(void *ptr, int histogram, int * array, int size)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            
            if (histogram >= ref->get_numberHisto_int())
                return -1;
            
            int * in_array = (ref->get_histo_array_int(histogram));
            if (array == NULL )
                return -1;
            else
            {
                // copy memory to avoid memory leaks...
                // numpy will handle its mem, we handle ours
                memcpy ( array, in_array , sizeof(int32_t)*size);
                
                delete [] in_array;
                return 0;
            }
        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }
    
    // get time resolution in ns per bin
    double get_binWidth_ns(void *ptr)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            return ref->get_binWidth_ns();
        }
        catch(...) // Catches everything
        {
           return -1.; //assuming -1 is an error condition. 
        }
    }

    // get t0 in ns 
    double get_t0_double(void *ptr, int histogram)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            return ref->get_t0_double(histogram);
        }
        catch(...) // Catches everything
        {
           return -1.; //assuming -1 is an error condition. 
        }
    }
    
    // get sample string
    int get_sample(void *ptr, char * sampleName)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            std::string res = ref->get_sample();
	    memcpy (sampleName, res.c_str(), sizeof(res.c_str()));
	    return 0;
        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }
    
    // get field string
    int get_field(void *ptr, char * fieldName)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            std::string res = ref->get_field();
	    memcpy (fieldName, res.c_str(), sizeof(res.c_str()));
	    return 0;
        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }

    // get orient string
    int get_orient(void *ptr, char * orientName)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            std::string res = ref->get_orient();
	    memcpy (orientName, res.c_str(), sizeof(res.c_str()));
	    return 0;
        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }

    // get temp string
    int get_temp(void *ptr, char * tempName)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            std::string res = ref->get_temp();
	    memcpy (tempName, res.c_str(), sizeof(res.c_str()));
	    return 0;
        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }

    // RDR attempt to add one new item  based on example get_temperatures_vector (changing double -> int)
    // get get_eventsHisto_vector
    int get_eventsHisto_vector(void *ptr, void * array)
// the logic is that this routine employs the previous 
// get_numberHisto_int to assess size 
// and produces the output as a pointer to the argument array, 
// suited for a numpy int array, 
// the calling python routines allocates it before the call
// (using the int for error management)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            
            int size = ref->get_numberHisto_int();
	    std::vector<long> in_array;  // top definition // capisco che in_array qui è definito come long, i.e. 64 bit 
            if (size==-1)
                return -1;
            else
                in_array = (ref->get_eventsHisto_vector()); // get vector from Alex if size is non zero

            if (array == NULL )
                return -1;
            else
            {
                // copy memory to avoid memory leaks...
                // numpy will handle its mem, we handle ours (vuol dire che il dtype di np definirà la dimensione delle words di python?)
                memcpy ( array, &in_array[0] , sizeof(in_array)*size); // a partire da &in_array[0] legge size long words e le passa in array
                // notice the way to copy from address of first vector element
                // purtroppo python si lamenta che ArgumentError: argument 2: <class 'TypeError'>: array must have data type int32
                return 0;
                // no need to delete memory allocation
                // because a vector is passed through the stack
            }
        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }

    // end RDR attempt

    // get comment string
    int get_comment(void *ptr, char * comment)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            std::string res = ref->get_comment();
	    memcpy (comment, res.c_str(), 62); // sizeof(res.c_str()));
	    return 0;
        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }

    // get run number int // another RDR addition
    int get_runNumber_int(void *ptr)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            return ref->get_runNumber_int();
        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }

    // get start time and date vector
    int get_timeStart_vector(void *ptr, char * timestring)
// produces the output as a pointer to the argument array, 
// suited for a string variable, 
// the calling python routines allocates it before the call
// (using the int for error management)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            
	    std::vector<std::basic_string<char> >  in_vector;  // top definition
            in_vector = ref->get_timeStart_vector();  // get vector from Alex
            std::string res = ""; 
            res += in_vector[0];  // convert to string time
            res += ' '; // spacer
            res += in_vector[1];  // convert to string date
            memcpy ( timestring, res.c_str() , sizeof(char)*18); // copy time+date
                                                            // to python string

            return 0;

        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }

    // get stop time and date vector
    int get_timeStop_vector(void *ptr, char * timestring)
// produces the output as a pointer to the argument array, 
// suited for a string variable, 
// the calling python routines allocates it before the call
// (using the int for error management)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            
	    std::vector<std::basic_string<char> >  in_vector;  // top definition
            in_vector = ref->get_timeStop_vector();  // get vector from Alex
            std::string res = ""; 
            res += in_vector[0];  // convert to string time
            res += ' '; // spacer
            res += in_vector[1];  // convert to string date
            memcpy ( timestring, res.c_str() , sizeof(char)*18); // copy time+date
                                                            // to python string

            return 0;

        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }

    // get length of temperatures vector
    int get_numberTemperature_int(void *ptr)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            return ref->get_numberTemperature_int();
        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }

    // get temperatures vector
    int get_temperatures_vector(void *ptr, double * array)
// the logic is that this routine employs the previous 
// get_numberTemperature_int to assess size 
// and produces the output as a pointer to the argument array, 
// suited for a numpy float64 array, 
// the calling python routines allocates it before the call
// (using the int for error management)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            
            int size = ref->get_numberTemperature_int();
	    std::vector<double> in_array;  // top definition
            if (size==-1)
                return -1;
            else
                in_array = (ref->get_temperatures_vector());

            if (array == NULL )
                return -1;
            else
            {
                // copy memory to avoid memory leaks...
                // numpy will handle its mem, we handle ours
                memcpy ( array, &in_array[0] , sizeof(double)*size);
                // notice the way to copy from address of first vector element
                return 0;
                // no need to delete memory allocation
                // because a vector is passed through the stack
            }
        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }

    // get std deviation of temperatures vector
    int get_devTemperatures_vector(void *ptr, double * array)
// the logic is that this routine employs the previous 
// get_numberTemperature_int to assess size 
// and produces the output as a pointer to the argument array, 
// suited for a numpy float64 array, 
// the calling python routines allocates it before the call
// (using the int for error management)
    {
        try
        {
            MuSR_td_PSI_bin * ref = reinterpret_cast<MuSR_td_PSI_bin *>(ptr);
            
            int size = ref->get_numberTemperature_int();
	    std::vector<double> in_array;  // top definition
            if (size==-1)
                return -1;
            else
                in_array = (ref->get_devTemperatures_vector());

            if (array == NULL )
                return -1;
            else
            {
                // copy memory to avoid memory leaks...
                // numpy will handle its mem, we handle ours
                memcpy ( array, &in_array[0] , sizeof(double)*size);
                // notice the way to copy from address of first vector element
                return 0;
                // no need to delete memory allocation
                // because a vector is passed through the stack
            }
        }
        catch(...) // Catches everything
        {
           return -1; //assuming -1 is an error condition. 
        }
    }

 
} //End C linkage scope.
