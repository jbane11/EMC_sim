{ char* macros[] =
    {   //list of scripts to load, end with '\0'
	  "av18_table.c",
	  //"sim_one.c",
      "\0"
    };
	
  if (*macros[0]!=0)
    cout << "\nrootlogon.C: Loading macro's:" << endl;
  for(UInt_t i=0; *macros[i]!=0; i++) {
    cout << "\t " << macros[i] << endl;
    gROOT->LoadMacro(macros[i]);
  }

  printf("\nrootlogon.C: Done!\n\n");
}
