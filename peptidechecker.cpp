// C++ program to check for side products of a peptide synthesis of a certain
//mass
//Author: Evan Munro (munro@stanford.edu)
//Date: July 29, 2019

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
using namespace std;

double sum_vector(vector<double>& masses) {
    double sum  = 0.0;
    for (int i = 0 ; i < masses.size(); i++){
        sum+= masses[i];
    }
    return sum;
}

//recursive function to search for side products of a peptide synthesis of
//a certain mass
void findCandidatePeptides(double target_mass,vector<char>& component_aminos,
                                vector<double>& component_masses, int n) {
    double water_mass = 18.0106;
    double cap_mass = 41.0265;

    double candidate_mass = sum_vector(component_masses);
    candidate_mass = candidate_mass - water_mass*(component_masses.size()-1) + cap_mass;

    if(abs(candidate_mass-target_mass)<1){
        cout << "Mass: "<< candidate_mass << " ";
        string s(component_aminos.begin(),component_aminos.end());
        cout << s << "\n";
    }

    if(candidate_mass < target_mass) {
        return;
    }

    for (int i = n ; i < component_masses.size(); i++) {
        vector<double> drop_one_mass(component_masses.size());
        copy(component_masses.begin(), component_masses.end(),drop_one_mass.begin());
        vector<char> drop_one_acids(component_aminos.size());
        copy(component_aminos.begin(), component_aminos.end(),drop_one_acids.begin());
        drop_one_mass.erase(drop_one_mass.begin()+i);
        drop_one_acids.erase(drop_one_acids.begin()+i);
        findCandidatePeptides(target_mass,drop_one_acids,drop_one_mass,i);
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
    int mass;
    //int min_mass;
    cout << "Enter the peptide string: ";
    cin >> peptide;
    vector<char> pept(peptide.begin(),peptide.end());
    vector<double> arr(pept.size());
    for (int i = 0; i < (pept.size()); i++) {
        arr[i] = amino_mass[pept[i]];
    }
    double total_mass = sum_vector(arr) - 18.0106*(arr.size()-1)+41.0265;
    cout << "Total Peptide Mass: " << total_mass << "\n";
    cout << "Enter the assumed mass of side product derived from LCMS output: ";
    cin >> mass;

    findCandidatePeptides(mass,pept,arr,0);

    return 0;
}
