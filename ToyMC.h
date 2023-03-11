#include <iostream>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <TH2F.h>
#include <TH2D.h>
#include <TFile.h>
#include <TString.h>
#include <TRandom3.h>
 #include <ctime>

using namespace std;

#define MAX_CTAG 4000

const double adt = 1.0;

class Hit {
  public:
    double time=-1;
    double energy=-1;    
    double start_time=-1;
    
    long seed=-1;
    int fill=-1;    
    int evt=-1;

    void Print(string v="") {        
        cout << v << "  time=" << time << "    energy=" << energy <<"     start time="<< start_time << "      fill=" << fill << "    random seed=" << seed<<endl;        
    }
};

// old clustering
bool Merge(Hit * this_cluster, Hit * this_hit, Hit * next_hit) {

    // this_cluster->Print();
    // this_hit->Print();    
    // next_hit->Print();

    bool merged = false;
    if(this_hit->time - this_cluster->start_time < adt) {    
        
        double merge_energy = this_cluster->energy + this_hit->energy;
        double merge_time = ( this_cluster->time*this_cluster->energy + this_hit->time*this_hit->energy ) / merge_energy;
        this_cluster->time = merge_time;
        this_cluster->energy = this_cluster->energy + this_hit->energy;
        merged = true;
    }
    else{
        merged = false;
    }

    return merged;

}

// time partioning

class ToyMC {
public:
    ToyMC(int ctag,int seed) {
        m_ctag = ctag;
        m_seed = seed;
        m_random = new TRandom3();        
    }

    void BookFileHists(TFile * file) {
        m_file = file;
        _BuildHistograms();        
    }

    void BuildRaw(int Nfills) {
        
        
        long start = time(0);
        for(m_this_fill=0;m_this_fill<Nfills;m_this_fill++) {
            m_this_seed = (m_ctag*1e3 + m_seed)*1e6+m_this_fill;
            m_random -> SetSeed(m_this_seed);

            if(m_this_fill%20000==0) {
                long now = time(0);
                cout << "processed "<< m_this_fill << " / " << Nfills << "\t\ttime="<< now-start << "s" << endl;
            }            
            _BuildTruth();
            _BuildRecon();
            // m_raw_hits[m_N_reco-1].Print("last hit:");
            _FillHistograms();
        }
    }

    void Finalize() {
        for(auto item : m_hists_1d){
            item.second -> Write();
        }
        for(auto item : m_hists_2d){
            item.second -> Write();
        }
        m_file->Close();
    };


private:
    int m_seed;
    
    long m_this_seed;
    int m_this_fill;

    int m_ctag;
    int m_N_reco;

    TFile * m_file;
    TRandom3 * m_random;

    Hit m_truth_hits[MAX_CTAG];
    Hit m_raw_hits[MAX_CTAG];

    map<string,TH1D*> m_hists_1d;
    map<string,TH2D*> m_hists_2d;

    void _BuildHistograms() {
        m_hists_2d.clear();
        m_hists_2d["raw"] = new TH2D(Form("Raw_ctag%d",m_ctag),"",500e3,0,500e3,10,0,10);
        m_hists_2d["truth"] = new TH2D(Form("Truth_ctag%d",m_ctag),"",500e3,0,500e3,10,0,10);

        m_hists_2d["raw_dt"] = new TH2D(Form("Raw_ctag%d_dt",m_ctag),"",500e3,225e3,275e3,100,0,10);        
        m_hists_2d["truth_dt"] = new TH2D(Form("Turth_ctag%d_dt",m_ctag),"",500e3,225e3,275e3,100,0,10);

        for(auto item : m_hists_1d){
            item.second -> SetDirectory(m_file);
        }
        for(auto item : m_hists_2d){
            item.second -> SetDirectory(m_file);
        }
    }

    void _BuildTruth() {
        // build truth
        for(int n=0;n<m_ctag;n++){             
            Hit hit;
            hit.time = m_random->Uniform(0,500e3);
            hit.start_time = hit.time;
            hit.seed = m_this_seed;
            hit.fill = m_this_fill;

            hit.energy = 1;
            m_truth_hits[n] = hit;
        }

        // sort truth by time
        std::sort(m_truth_hits, m_truth_hits + m_ctag,
          [](Hit const & a, Hit const & b) -> bool
          { return a.time < b.time; } );
    }

    void _BuildRecon() {
        // build recon from truth
        Hit this_recon;
        int n_reco = 0;
        int n_truth = 0;
        while(n_truth<m_ctag) {
            // fill this_recon
            if(this_recon.time<0) {
                this_recon = m_truth_hits[n_truth];                
            }
            // compare this_recon to this_hit and next_hit
            else {
                bool merged = Merge(&this_recon, m_truth_hits+n_truth, m_truth_hits+n_truth+1);
                if(!merged){
                    //generate a raw hit
                    m_raw_hits[n_reco++] = this_recon;
                    this_recon = m_truth_hits[n_truth];
                }
            }
            n_truth++;
        }
        //fill last hit
        m_raw_hits[n_reco++] = this_recon;
        m_N_reco = n_reco;
    }

    void _FillHistograms() {
        // reco
        for(int n=0;n<m_N_reco;n++) {
            m_hists_2d["raw"] -> Fill(m_raw_hits[n].time,m_raw_hits[n].energy);
            if(n>0) {
                double dt = m_raw_hits[n].time-m_raw_hits[n-1].time;
                m_hists_2d["raw_dt"] -> Fill(m_raw_hits[n].time,dt);
            }
        }

        //truth
        for(int n=0;n<m_ctag;n++) {
            m_hists_2d["truth"] -> Fill(m_truth_hits[n].time,m_truth_hits[n].energy);
            if(n>1) {                
                m_hists_2d["truth_dt"] -> Fill(m_truth_hits[n].time,m_truth_hits[n].time-m_truth_hits[n-1].time);
            }
        }
    }
};