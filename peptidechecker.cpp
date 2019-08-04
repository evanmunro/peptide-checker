// C++ program to check for side products of a peptide synthesis of a certain
//mass
//Author: Evan Munro (munro@stanford.edu)
//Date: July 29, 2019

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
using namespace std;

double CAP_MASS = 41.0265;
double WATER_MASS = 18.0106;
double sum_vector(vector<double>& masses) {
    double sum  = 0.0;
    for (int i = 0 ; i < masses.size(); i++){
        sum+= masses[i];
    }
    return sum;
}

//recursive function to search for side products of a peptide synthesis of
//a certain mass
void findCandidatePeptides(vector<char>& orig_pep,double target_mass,vector<char>& component_aminos,
                                vector<double>& component_masses, int end_ignore, int n) {
    double tolerance = 1.0;
    double candidate_mass = sum_vector(component_masses);
    candidate_mass = candidate_mass - WATER_MASS*(component_masses.size()-1) + CAP_MASS;

    if(abs(candidate_mass-target_mass)<tolerance){
        cout << candidate_mass << ", ";
        vector<char> side_product(orig_pep.begin(),orig_pep.end());
        int j = orig_pep.size()-1;
        int i = component_aminos.size()-1;
        while(j>=0 && i >=0 ) {
            if(component_aminos[i]== orig_pep[j]){
                j= j-1;
                i=i-1;
            }
            else{
                side_product[j] = '*';
                j= j-1;
            }
        }
        while(j>=0) {
            side_product[j] = '*';
            j= j-1;
        }
        string s(side_product.begin(),side_product.end());
        cout << s << "\n";
    }

    if(candidate_mass < target_mass) {
        return;
    }
    for (int i = n ; i < component_masses.size()-end_ignore; i++) {
        vector<double> drop_one_mass(component_masses.size());
        copy(component_masses.begin(), component_masses.end(),drop_one_mass.begin());
        vector<char> drop_one_acids(component_aminos.size());
        copy(component_aminos.begin(), component_aminos.end(),drop_one_acids.begin());
        drop_one_mass.erase(drop_one_mass.begin()+i);
        drop_one_acids.erase(drop_one_acids.begin()+i);
        findCandidatePeptides(orig_pep,target_mass,drop_one_acids,drop_one_mass,end_ignore,i);
    }



}


int main()
{
    map<char, double> amino_mass = {
        { 'A', 89.0477 },
        { 'R', 174.1117 },
        { 'N', 132.0535},
        {'D', 133.0375},
        { 'C', 121.0197},
        { 'E', 147.0532 },
        { 'Q', 146.0691 },
        {'G', 75.032},
        { 'H', 155.0695 },
        { 'I', 131.0946 },
        { 'L', 131.0946},
        {'K', 146.1055},
        { 'M', 149.051 },
        { 'F', 165.079 },
        { 'P', 115.0633},
        {'S', 105.0426},
        { 'T', 119.0582},
        { 'W', 204.0899 },
        { 'Y', 181.0739 },
        {'V', 117.079}
    };
    string peptide;
    double mass;
    int end_ignore;
    //int min_mass;
    cout << "Enter the peptide string: ";
    cin >> peptide;
    vector<char> pept(peptide.begin(),peptide.end());
    vector<double> arr(pept.size());
    for (int i = 0; i < (pept.size()); i++) {
        arr[i] = amino_mass[pept[i]];
    }
    double total_mass = sum_vector(arr) - WATER_MASS*(arr.size()-1)+CAP_MASS;
    cout << "Total Peptide Mass: " << total_mass << "\n";
    cout << "Enter the assumed mass of unknown side product derived from LCMS output,"
    " \n or 0 if you want to just see a list of masses for truncated peptides: ";
    cin >> mass;
    //Assuming synthesis starting from the end of the peptide string
    if (mass ==0.0) {
        double running_sum = arr[pept.size()-1];
        string running_pept= peptide.substr(pept.size()-1,1);
        for (int i =pept.size()-2; i >=0 ; i=i-1) {
            running_sum += arr[i]-WATER_MASS;
            running_pept = peptide[i] + running_pept;
            cout << running_sum+CAP_MASS << " " << running_pept << "\n";
        }
        return 0;
    }
    else {
        cout << "Enter the number of amino acids at the end of"
        " the peptide, \n that you would like to assume do not fail to couple when checking combinations: ";
        cin >> end_ignore;
        findCandidatePeptides(pept,mass,pept,arr,end_ignore,0);
    }
    return 0;
}
