type file;

app (file oSG, file eSG) runSpotsGen (string fltfldr, string ftm, string paramsfn, file DummyStepD[])
{
 SpotsGen fltfldr ftm paramsfn stdout=@filename(oSG) stderr=@filename(eSG);
}

app (file oGPSO, file eGPSO) runGenPSOut (string foldr, string filestm, string tempofldr, int iternr, string parafn, int RingNr, file DummyStepB)
{
 GenPSOut foldr filestm tempofldr iternr parafn RingNr stdout=@filename(oGPSO) stderr=@filename(eGPSO);
}

app (file oCR, file eCR) runCalcRad (string fldr, string fstm, string parfn, int Ring, file DummyStepC[])
{
 CalcRad fldr fstm parfn Ring stdout=@filename(oCR) stderr=@filename(eCR);
}

app (file oPPF, file ePPF) runPeaksPerFile (string foldern, string fs, int fnr, string tmpfldr, string drk, string params, int Rin)
{
 PeaksPerFile foldern fs fnr tmpfldr drk params Rin stdout=@filename(oPPF) stderr=@filename(ePPF);
}

app (file oFO, file eFO) runFindOverlap (string fldrn, string fstem, string tmpfl, string paramsf, int RNr, file DummystepA[])
{
 FindOverlap fldrn fstem tmpfl paramsf RNr stdout=@filename(oFO) stderr=@filename(eFO);
}

# Parameters to be modified #############

string folder = @arg("folder","/data/tomo1/Sharma_Jun13/BINFiles/");
int startnr = @toInt(@arg("startnr","1"));
int endnr = @toInt(@arg("endnr","600"));
int Rings[] = [1, 2, 3, 4];
string darkfn = @arg("dark","/data/tomo1/Sharma_Jun13/BINFiles/Dark_000001.bin");
string parameterfilename = @arg("paramsfile","/clhome/TOMO1/PeaksAnalysisHemant/PeaksFittingCode/ParamsFile.txt");
string filestem = @arg("filestem","93_44_FF_Volume_Scan_Final_10");
string tempfolder = @arg("tempfolder","Temp2");

# End parameters ########################

string PreFix0 = @strcat("CalcRad_",filestem);
file simDout[]<simple_mapper;location="output",prefix=PreFix0,suffix=".out">;
file simDerr[]<simple_mapper;location="output",prefix=PreFix0,suffix=".err">;
foreach RingsNrs in Rings {
    # Find peaks per file ###################
    string PreFix1 = @strcat("PeaksPerFile_",RingsNrs,"_");
    file simAerr[]<simple_mapper;location="output",prefix=PreFix1,suffix=".err">;
    file simAout[]<simple_mapper;location="output",prefix=PreFix1,suffix=".out">;
    foreach i in [startnr:endnr] {
       (simAout[i],simAerr[i]) = runPeaksPerFile(folder,filestem,i,tempfolder,darkfn,parameterfilename,RingsNrs);
    }
    # Find overlapping files ################
    file simBout<single_file_mapper; file=@strcat("output/FindOverlap_",filestem,RingsNrs,".out")>;
    file simBerr<single_file_mapper; file=@strcat("output/FindOverlap_",filestem,RingsNrs,".err")>;
    (simBout,simBerr) = runFindOverlap(folder,filestem,tempfolder,parameterfilename,RingsNrs,simAerr);
    # Generate PS Output ####################
    string PreFix2 = @strcat("GeneratePSOut_",RingsNrs,"_");
    file simCerr[]<simple_mapper;location="output",prefix=PreFix2,suffix=".err">;
    file simCout[]<simple_mapper;location="output",prefix=PreFix2,suffix=".out">;
    int NumberOfFiles=(endnr - startnr + 1);
    foreach i in [1:NumberOfFiles] {
        (simCout[i],simCerr[i]) = runGenPSOut(folder,filestem,tempfolder,i,parameterfilename,RingsNrs,simBout);
    }
    # Calculate Radii #######################
    (simDout[RingsNrs],simDerr[RingsNrs]) = runCalcRad(folder,filestem,parameterfilename,RingsNrs,simCerr);
}

# Generate input files ##################
file simEout<single_file_mapper; file=@strcat("output/TiltBCLsd_",filestem,".out")>;
file simEerr<single_file_mapper; file=@strcat("output/TiltBCLsd_",filestem,".err")>;
(simEout,simEerr) = runSpotsGen(folder,filestem,parameterfilename,simDout);
