
#include "Trait.h"

using namespace MeshLib;


void Trait::clear( Trait * & pT )
{
    if ( pT == NULL ) return;

    clear( pT->next() );
    delete( pT );
    pT = NULL;  
};

/*! Modification on July 26, 2006.  To discriminate the partial string case.
e.g. sfather=( and r=( */
void Trait::updateTraitString(std::string &traitString, std::string &traitName, std::string &traitValue)
{
    int sp, ep, sfp;
    sfp = (int)traitString.find(" "+traitName+"=(");
    sp = (int)traitString.find(traitName+"=(");
    if ( sp==0 || sfp != (int)std::string::npos )
    {
        int temp;
        if ( sp == 0 ) temp = sp;
        else temp = sfp;

        sp = (int)traitString.find_first_of("(", temp);
        ep = (int)traitString.find_first_of(")", temp);
        if ( sp != (int)std::string::npos && ep != (int)std::string::npos )
        {
            // get inside of the parenthesis
            sp++;
            ep--;
            traitString.replace(sp, ep-sp+1, traitValue);
        }
        else
        {
            // string is corrupt
            std::cerr << "ERROR: Trait string missing parenthesis: " << traitString << std::endl;
            return;
        }
    }
    else
    {
        // trait not present in string
        if ( traitString.length() )
        {
            traitString = traitName + "=(" + traitValue + ") " + traitString;
        }
        else
        {
            traitString = traitName + "=(" + traitValue + ")";
        }

        return;
    }
}

// see header for comments
std::string Trait::getTraitValue(std::string &traitString, std::string &traitName)
{
    int sp, ep, sfp;
    sfp = (int)traitString.find(" "+traitName+"=(");
    sp = (int)traitString.find(traitName+"=(");
    if ( sp == 0 || sfp != (int)std::string::npos )
    {
        int temp;
        if ( sp == 0 ) temp = sp;
        else temp = sfp;

        std::string isolatedTraitString(traitString.substr(temp));
        sp = (int)isolatedTraitString.find_first_of("(");
        ep = (int)isolatedTraitString.find_first_of(")");
        if ( sp != (int)std::string::npos )
        {
            // the value was found
            sp++;
            return isolatedTraitString.substr(sp, ep-sp);
        }
        else
        {
            std::cerr << "ERROR: String missing parenthesis: " << isolatedTraitString << std::endl;
            return std::string();
        }
    }

    return std::string();
}
