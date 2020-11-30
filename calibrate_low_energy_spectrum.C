//Program to fit peaks in low energy 25Mg(d,p)26Mg data, fit a polynomial and output fit parameters


void calibrate_low_energy_spectrum(){
  cout << "HEllow orld" << endl;

  int run_array[16] = {6, 9, 11, 14, 17, 23, 25, 26, 27, 30, 31, 34, 96, 99, 102, 104};
  int angle_array[16] = {10, 7, 13, 16, 29, 30, 33, 36, 39, 22, 25, 28, 45, 48, 52, 55};
  TFile *f = new TFile("all_spectra.root");

  int no_calibration_peaks = 8;

  double xbin, ybin, xmin, xmax, ymin, ymax, include_value;
  xbin=4096; ybin=100; xmin=0; xmax=4095; ymin=0; ymax=0;


  ifstream calibration_peaks_file;
  string calibration_peaks_string;
  double start, end, ex_energy;
  TH1F *h1 = new TH1F("h1","h1 title",xbin,xmin,xmax);

  int no_peaks = 8;
  double centroid_array[no_peaks];
  double ex_energy_array[no_peaks];
  double centroid_uncert_array[no_peaks];
  double ex_energy_uncert_array[no_peaks];


  int temp_index = 6;

  for(int i = temp_index; i < (temp_index + 1); i++){ // Iterating over runs
    cout << run_array[i] << " " << angle_array[i] << endl;
    h1 = (TH1F*)f->Get(Form("run%i", run_array[i]));
    h1->Draw();

    //Open file containing start parameter, end parameter and excitation energy
    calibration_peaks_file.open(Form("low_energy_calibration_peaks/run_%i_calibration.txt", run_array[i]));
    cout << "opening file " << Form("low_energy_calibration_peaks/run_%i_calibration.txt", run_array[i]) << endl;

    for(int j = 0; j < no_peaks; j++){ //Iterating over peaks in a run
      cout << j << endl;
      getline(calibration_peaks_file, calibration_peaks_string);
      const char *calib_char = calibration_peaks_string.c_str();
      cout << calib_char << endl;
      sscanf(calib_char, "%lf\t%lf\t%lf", &ex_energy, &start, &end);
      cout << start << " " << end << " " << ex_energy << endl;



      TF1 *f1;
      f1 = new TF1("f1", "gaus", start, end);
      h1->Fit("f1", "R", "", start, end);
      centroid_array[j] = f1->GetParameter(1);
      centroid_uncert_array[j] = f1->GetParError(1);

      cout << f1->GetParameter(1) << " " <<  centroid_array[j] << endl;
      ex_energy_array[j] = ex_energy;
      cout << ex_energy << " " << ex_energy_array[j] << endl;

      ex_energy_uncert_array[j] = 0.0001;

    }


    // TGraph *gr1 = new TGraph(no_peaks, centroid_array, ex_energy_array);
    // gr1->Draw("A*");

    TGraphErrors *gr1_errors = new TGraphErrors(no_peaks, centroid_array, ex_energy_array, centroid_uncert_array, ex_energy_uncert_array);
    gr1_errors->Draw("A*");


    TF1 *channel_to_exc = new TF1("channel_to_exc", "pol 2", 500, 4000);
    channel_to_exc->SetParameter(0, 7.59);
    channel_to_exc->SetParameter(1, -1.01971e-03);
    channel_to_exc->SetParameter(2, -3.55763e-08);

    // gr1->Fit(channel_to_exc, "R");


    gr1_errors->Fit(channel_to_exc, "R");
    cout << "Chi-sqaured is: " << channel_to_exc->GetChisquare() << endl;





    calibration_peaks_file.close();

    //Outputting fit parameters to file


    ofstream channel_to_ex_fit_parameters_output;
    // channel_to_ex_fit_parameters_output.open(Form("exc_to_channel_parameters/exc_to_channel_parameters_run%i.txt", run_array[i]));
    channel_to_ex_fit_parameters_output.open(Form("exc_to_channel_parameters/exc_to_channel_parameters_angle%i.txt", angle_array[i]));

    channel_to_ex_fit_parameters_output << channel_to_exc->GetParameter(0) << "\t"
    << channel_to_exc->GetParameter(1) << "\t"
    << channel_to_exc->GetParameter(2) << endl;

    channel_to_ex_fit_parameters_output << channel_to_exc->GetParError(0) << "\t"
    << channel_to_exc->GetParError(1) << "\t"
    << channel_to_exc->GetParError(2) << endl;


    cout << "Outputted to " << Form("exc_to_channel_parameters/exc_to_channel_parameters_angle%i.txt", angle_array[i]) << endl;

    channel_to_ex_fit_parameters_output.close();



  }

  // TH1F *h1 = new TH1F("h1","h1 title", xbin,xmin,xmax); //setting up histogram to extract from root file
  // h1 = (TH1F*)f->Get(Form("run%d", calib_run_no));


}
