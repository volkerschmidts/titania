**tnafiledata[*'FilePath'*]**: full file path  
**tnafiledata[*'Title'*]**: trimed and split input name
**tnafiledata[*'NOR']**: number of RDCs  
**tnafiledata[*'NOA'*]**: number of Atoms  
**tnafiledata[*'NOS'*]**: number of RDC sets  
**tnafiledata[*'Sets'*]**:
- [string] RDC sets  

**tnafiledata[*'RDCs'*]**:  
- [string] RDCs  

**tnafiledata[*'SECONDA'*]**:  
- dict  
    - Cutoff
    - Eigenmodes[]   
    - SECONDA[]   
        - Index  
        - EigVal  
        - Kappa  
        - csum  

**tnafiledata[*'Het_csum'*]**:  
- [double] a^2  

**tnafiledata[*'Qfactor'*]**:  
- [Set][Iteration] double  

**tnafiledata[*'ChiralVolume'*]**:  
- dict
    - Centers[]
    - Volumes[Centers] float  

**tnafiledata[*'PolarAngles'*]**:  
- dict  
    - RDC(dict)  
        - inital(dict): theta, phi  
        - final (dict): theta, phi  

**tnafiledata[*'Dynamics'*]**:  
- dict  
    - Soverall  
    - Center(dict)  
        - Srdc  
        - Sax  
        - eta  
        - phi  
        - chi2  
        - Dmax  

**tnafiledata[*'Iteration'*]**: Saves current Iteration  


**tnafiledata[*'Redundants'*]**:  
- dict  
    - Iteration[]  
    - Delta\_x[]  
    - Delta\_p[]  
    - StopCrit[]  
    - Damping(dict)  
        - effectiv(dict)  
            - Overall[]  
            - Bond[]   
            - Angle[]  
            - Dihedral[]  
            - RDC[]  
            - ChiralVolume[]  
        - static(dict)  
            - Bond  
            - Angle  
            - Dihedral  
            - RDC  
            - ChiralVolume  

**tnafiledata[*'StopCrit'*]**:  
- dict
    - Values(dict)
        - SaupeMean[]  
        - SaupeSigma[]  
        - PolarMean[]   
        - PolarSigma[]  
        - VectorLength[]  
        - Qfactor[]  
    - Thresholds(dict)
        - SaupeMean  
        - SaupeSigma  
        - PolarMean   
        - PolarSigma  
        - VectorLength  
        - Qfactor  


