//Macro to simply get histogram and zoom in on a small area of it
//has been edited to fit 5.7 MeV doublet and calculate excitation energy of peaks
const char* ConvertDoubleToString(double value){ std::stringstream ss ;
  ss << value;
  const char* str = ss.str().c_str();
  return str;
}



void zoom(){
  #include <vector>

  //Booleans
  bool doublet_boolean = 1;
  bool batch_boolean = 1;

  if(batch_boolean){
    gROOT->SetBatch(kTRUE);
  }
  else{
    gROOT->SetBatch(kFALSE);
  }


  cout <<" Hello" << endl;


  //Defining variables
  double xbin, ybin, xmin, xmax, ymin, ymax, include_value;
  xbin=4096; ybin=100; xmin=0; xmax=4095; ymin=0; ymax=0;

  int start, end, run_file_counter, corresp_angle, file_flag, i;
  int run_no=0;
  int total_runs=107;
  string s;

  std::vector<int> corresp_angle_vector;
  std::vector<double> corresp_angle_double_vector;
  std::vector<double> centroid_vector;
  std::vector<double> centroid_diff_vector;

  // std::TVectorD<int> corresp_angle_vector;
  // std::vector<double> centroid_diff_vector;

  TFile *f = new TFile("all_spectra.root");


  TH1F *h1 = new TH1F("h1","h1 title",xbin,xmin,xmax);

  ifstream peaklocationfile;

  //Getting histogram

  ofstream outputfile;
  outputfile.open("zoom_outputs.txt");

  ofstream exc_calib_output_file;
  exc_calib_output_file.open("exc_calib_output.txt");

  double doublet_large_peak_centroid, doublet_large_peak_height, doublet_large_peak_width;
  double doublet_small_peak_centroid, doublet_small_peak_height, doublet_small_peak_width;

  cout << "Conor 1" << endl;
  //Run loop begins
  for(run_no=0; run_no < 105; run_no++){
    //peaklocationfile.open("4+_4.901_peaklocations.txt");
    peaklocationfile.open("5.7_total_peaklocations.txt");
    // peaklocationfile.open("3+_6.125_peaklocations.txt");
    // peaklocationfile.open("2+_5.292_peaklocations.txt");
    // peaklocationfile.open("4+_5.474_peaklocations.txt");


    cout << "Conor 2" << endl;

    // printf("run number %d\n",run_no);
    h1 = (TH1F*)f->Get(Form("run%d", run_no));

    file_flag=0;
    //Getting details from input file
    for (i = 0; i < total_runs; i++){
      // printf("loop number %d\n",i);

      getline(peaklocationfile,s);
      const char *s_char = s.c_str();
      corresp_angle=0;
      include_value=0;
      // cout << s_char << endl;
      if(doublet_boolean){
        sscanf(s_char,"%d\t%d\t%d\t%d\t%lf",&run_file_counter,&start,&end, &corresp_angle, &include_value); //picking out variables
      }
      else{
        sscanf(s_char,"%d %d %d %d %lf",&run_file_counter,&start,&end, &corresp_angle, &include_value); //picking out variables
      }

      // cout << run_file_counter << run_no << endl;
      if(run_file_counter == run_no){
        cout << "location of " << start << " " << end << " found for run " << run_no << endl;
        cout << corresp_angle << " " << include_value << endl;
      //   cout << "run " << run_no << " is a match, zooming in on it" << endl;
      //   file_flag=1;
      if(include_value==1){
        TCanvas *c1 = new TCanvas();
        h1->Rebin();
        // h1->Rebin();
        h1->GetXaxis()->SetRangeUser((start-150),(end+150));
        h1->Draw();
        TLine *start_line = new TLine(start,0,start,1000);
        TLine *end_line = new TLine(end,0,end,1000);
        start_line->SetLineColor(kRed);
        end_line->SetLineColor(kRed);
        start_line->Draw();
        end_line->Draw();

        //Fitting single gaus to get initial centroid value
        double large_peak_centroid, large_peak_height, large_peak_width;
        double centroid_diff, centroid_diff_uncert, centroid_diff_energy, centroid_diff_energy_uncert;
        TF1 *f1;
        f1 = new TF1("f1", "gaus", start, end);
        char pdf_name[20];
        h1->Fit("f1", "R", "", start, end);

        large_peak_height = f1->GetParameter(0);
        large_peak_centroid = f1->GetParameter(1);
        large_peak_width = f1->GetParameter(2);



        if(doublet_boolean){
          TPaveText *pt_single = new TPaveText(start + 100, f1->GetParameter(0)*0.5, end + 100, f1->GetParameter(0));
          char corresp_angle_char[40];
          sprintf(corresp_angle_char, "Lab angle: %i degrees", corresp_angle);
          pt_single->AddText(corresp_angle_char);
          pt_single->Draw();

          sprintf(pdf_name, "5.7_peak_fits/run%i_single_gaus_5.7.pdf", run_no);
          c1->SaveAs(pdf_name);
          delete pt_single;


          char double_gaus[16]; //Formula for double gaussian fit
          sprintf(double_gaus, "gaus(0)+gaus(3)");
          TF1* f1_double;
          f1_double = new TF1("f1_double", double_gaus, start, end);
          f1_double->SetParameters(large_peak_height, large_peak_centroid, 8,
            0.2 * large_peak_height, large_peak_centroid + 20, 8);
            f1_double->FixParameter(2, 6.5);
            f1_double->FixParameter(5, 6.5);

            printf("%lf %lf %lf\n", large_peak_height, large_peak_centroid, large_peak_width);
            double *double_params = f1_double->GetParameters();
            cout << "double gaus params: " << endl;
            for (int j = 0; j < 6; j++){
              cout << double_params[j] << endl;
            }
            h1->Fit("f1_double", "BR", "", start, end);

            TPaveText *pt_double = new TPaveText(start + 100, f1_double->GetParameter(0)*0.5, end + 100, f1_double->GetParameter(0));

            centroid_diff = f1_double->GetParameter(4) - f1_double->GetParameter(1);
            centroid_diff_uncert = sqrt(f1_double->GetParError(4) * f1_double->GetParError(4) + f1_double->GetParError(1) * f1_double->GetParError(1));


            doublet_large_peak_height = f1_double->GetParameter(0);
            doublet_large_peak_centroid = f1_double->GetParameter(1);
            doublet_large_peak_width = f1_double->GetParameter(2);

            doublet_small_peak_height = f1_double->GetParameter(3);
            doublet_small_peak_centroid = f1_double->GetParameter(4);
            doublet_small_peak_width = f1_double->GetParameter(5);



            pt_double->AddText(corresp_angle_char);

            //Adding centroid difference on to histogram
            // TLatex *centroid_diff_latex = new TLatex;
            // centroid_diff_latex->SetTextSize(100000);
            // centroid_diff_latex->SetTextAlign(12);
            // centroid_diff_latex->DrawLatex(end, 500, "test latex");
            // h1->Draw();
            // c1->Draw();

            char centroid_diff_char[40];
            sprintf(centroid_diff_char, "Centroid diff: \n%.1lf", centroid_diff);
            printf("Centroid diff is %lf\n", (centroid_diff));
            cout << centroid_diff_char << endl;
            pt_double->AddText(centroid_diff_char);


            pt_double->Draw();



            sprintf(pdf_name, "5.7_peak_fits/run%i_double_gaus_5.7.pdf", run_no);
            c1->SaveAs(pdf_name);

            corresp_angle_vector.push_back(corresp_angle);
            centroid_vector.push_back(large_peak_centroid);
            centroid_diff_vector.push_back(centroid_diff);
            corresp_angle_vector.push_back(double(corresp_angle));

            cout << corresp_angle_vector[0] << " " << centroid_diff_vector[0] << endl;
            outputfile << corresp_angle << " " << large_peak_centroid << " " << centroid_diff << " ";
        } //End of doublet boolean


        //Converting large peak centroid to excitation energy

        string channel_to_ex_fit_parameters_string, channel_to_ex_fit_parameters_uncert_string;
        double channel_to_ex_fit_parameter_a, channel_to_ex_fit_parameter_b, channel_to_ex_fit_parameter_c, exc_energy;
        double channel_to_ex_fit_parameter_a_uncert, channel_to_ex_fit_parameter_b_uncert, channel_to_ex_fit_parameter_c_uncert;


        //Opening up channel to exc calibration parameters file
        ifstream channel_to_ex_fit_parameters_input;
        channel_to_ex_fit_parameters_input.open(Form("exc_to_channel_parameters/exc_to_channel_parameters_angle%i.txt", corresp_angle));
        getline(channel_to_ex_fit_parameters_input, channel_to_ex_fit_parameters_string);
        const char* channel_to_ex_fit_parameters_char = channel_to_ex_fit_parameters_string.c_str();
        sscanf(channel_to_ex_fit_parameters_char, "%lf\t%lf\t%lf", &channel_to_ex_fit_parameter_a, &channel_to_ex_fit_parameter_b, &channel_to_ex_fit_parameter_c);

        getline(channel_to_ex_fit_parameters_input, channel_to_ex_fit_parameters_uncert_string);
        const char* channel_to_ex_fit_parameters_uncert_char = channel_to_ex_fit_parameters_uncert_string.c_str();
        sscanf(channel_to_ex_fit_parameters_uncert_char, "%lf\t%lf\t%lf", &channel_to_ex_fit_parameter_a_uncert, &channel_to_ex_fit_parameter_b_uncert, &channel_to_ex_fit_parameter_c_uncert);

        // sscanf(s_char,"%d\t%d\t%d\t%d\t%lf",&run_file_counter,&start,&end, &corresp_angle, &include_value); //picking out variables

        printf("Parameters: %lf\t%lf\t%lf\n", channel_to_ex_fit_parameter_a, channel_to_ex_fit_parameter_b, channel_to_ex_fit_parameter_c);
        printf("Parameters uncertainties: %lf\t%lf\t%lf\n", channel_to_ex_fit_parameter_a_uncert, channel_to_ex_fit_parameter_b_uncert, channel_to_ex_fit_parameter_c_uncert);



        channel_to_ex_fit_parameters_input.close();

        //For 5.7 MeV peak from single fit
        exc_energy = channel_to_ex_fit_parameter_a +
        channel_to_ex_fit_parameter_b * large_peak_centroid +
        channel_to_ex_fit_parameter_c * large_peak_centroid * large_peak_centroid;


        if(doublet_boolean){
          //For 5.7 MeV peak from doublet fit
          exc_energy = channel_to_ex_fit_parameter_a +
          channel_to_ex_fit_parameter_b * doublet_large_peak_centroid +
          channel_to_ex_fit_parameter_c * doublet_large_peak_centroid * doublet_large_peak_centroid;

          //For 5.69 MeV peak
          exc_energy = channel_to_ex_fit_parameter_a +
          channel_to_ex_fit_parameter_b * doublet_small_peak_centroid +
          channel_to_ex_fit_parameter_c * doublet_small_peak_centroid * doublet_small_peak_centroid;

          centroid_diff_energy = centroid_diff * -channel_to_ex_fit_parameter_b;
          // centroid_diff_energy_uncert = centroid_diff * channel_to_ex_fit_parameter_b_uncert;

          centroid_diff_energy_uncert = centroid_diff_energy * sqrt( (channel_to_ex_fit_parameter_b_uncert/channel_to_ex_fit_parameter_b) * (channel_to_ex_fit_parameter_b_uncert/channel_to_ex_fit_parameter_b) + (centroid_diff_uncert/centroid_diff) * (centroid_diff_uncert/centroid_diff) );


          outputfile << 1000.0 * centroid_diff_energy << " " << 1000.0 * centroid_diff_energy_uncert << endl;



        }

        exc_calib_output_file  << corresp_angle << "\t" << exc_energy << endl;



      } //End of include value boolean
        // break;
      } //End of run_file_counter == run_no boolean
    } //End of finding correct run no iterator
    peaklocationfile.close();
  } //End of run no iterator

  double* centroid_diff_array = &centroid_diff_vector[0];
  int* corresp_angle_array = &corresp_angle_vector[0];
  double* corresp_angle_double_array = &corresp_angle_double_vector[0];




  // for (int i = 0; i < 40; i++){
  //   cout << corresp_angle_array[i] <<" " << centroid_diff_array[i] << "\n";
  //   outputfile <<corresp_angle_array[i] << "\t" << centroid_diff_array[i] << "\n";
  // }
  outputfile.close();
  exc_calib_output_file.close();
  return;

  int plotting_points_no;
  // cout << corresp_angle_array.size() <<  endl;
  cout << corresp_angle_array[0] <<" " << typeid(*corresp_angle_array).name() << endl;
  // for (int i = 0; i < plotting_points_no; i++){
  //   cout <<"t"<< endl;
  // }
  plotting_points_no = 0;
  cout << "Output of begin and end: ";
  for (auto i = corresp_angle_vector.begin(); i != corresp_angle_vector.end(); ++i)
      cout << *i << "\n";
      cout << i << endl;
      plotting_points_no++;

  cout << "plotting points number: " << plotting_points_no << endl;
  // double corresp_angle_double_array[plotting_points_no];
  // int iterator_counter = 0;
  // for (auto i = corresp_angle_vector.begin(); i != corresp_angle_vector.end(); ++i)
  //     cout << *i << "\n";
  //     // corresp_angle_double_array[iterator_counter] = *i;
  //     cout << i << endl;
  //     cout <<corresp_angle_double_array[iterator_counter]<<endl;
  //     iterator_counter++;

  // cout << corresp_angle_double_array[0] << endl;
  //Plotting out centroid diff change across angles
  // TGraph *gr1 = new TGraph(corresp_angle_vector, centroid_diff_vector);
  TGraph *gr1 = new TGraph(plotting_points_no, corresp_angle_double_array, centroid_diff_array);
  gr1->Draw("A*");


  //Now rereading in output file and calculting mean and standard deviation of energy differences
  printf("Now rereading in output file and calculating mean and standard deviation of energy differences");

  int angle_array[16] = {7, 10, 13, 16, 29, 30, 33, 36, 39, 22, 35, 28, 45, 48, 52, 55};

  ifstream input_file;
  input_file.open("zoom_outputs.txt");

  if(doublet_boolean){
    for(int i = 0; i < 16; i++){
      cout << "angle is " << angle_array[i] << endl;

      for(int j = 0; j < 27; j++){

      }


    }

    input_file.close();



  }





  return;
}
