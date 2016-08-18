package gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import java.util.ArrayList;

public class GWASFormatter2 {

	/**
	 * 
	 * Attributes
	 */
	//private static HashMap <Integer, Integer> patients = new HashMap ();//<patientID, >
	@SuppressWarnings({ "rawtypes", "unchecked" })
	protected static HashMap <String, HashMap <Integer, String>> patient_snps = new HashMap ();//<snp, {PatientID:genotype, PatientID:genotype, ....}>
	@SuppressWarnings({ "unchecked", "rawtypes" })
	protected static HashSet <Integer> all_pt= new HashSet ();
	//protected static HashMap <Integer, ArrayList <String>> patient_snp_list= new HashMap ();//<PatientID, {SNP1, SNP2...}>
	@SuppressWarnings({ "unchecked", "rawtypes" })
	protected static HashMap <String, String> snp_refs = new HashMap ();//<SNP, Ref_genotype>
	@SuppressWarnings({ "unchecked", "rawtypes" })
	protected static HashMap <Integer, ArrayList <String>> data = new HashMap ();// <Patient, {GENOTYPE.....}, this comprises the TOTAL data
	@SuppressWarnings({ "unchecked", "rawtypes" })
	protected static HashMap <Integer, ArrayList <String>> metadata = new HashMap ();//<patient, {fam.id, ind.id, pat.id, mat.id, sex, pheno}>
	//mat and pat ids are defaulted at 0, phenotype is fusion factor and famID is patient id 
	@SuppressWarnings({ "unchecked", "rawtypes" })
	protected static HashMap <Integer, ArrayList <String>> md = new HashMap ();//patient --> {sex, phenotype}
	@SuppressWarnings({ "unchecked", "rawtypes" })
	protected static HashSet <String> unwanted = new HashSet ();
	@SuppressWarnings({ "unchecked", "rawtypes" })
	protected static ArrayList <String> all_snps = new ArrayList ();
	//protected static HashMap <String, Integer> snp_pos = new HashMap ();//<SNP, BP Position>
	protected static int total;
	protected static int size=patient_snps.size();
	protected static boolean header=false;

	//=+=+=+=+=+=+=+=+=+=+Usage=+=+=+=+=+=+=+=+=+=+=+=+=+=+
	public static void help (){
		System.out.println("Usage of GWASFormatter\n"
				+
				"Use java GWASFormatter2 <varScan file> <patient phenotype file> -hd/-nd threshold numThreads\n"
				+ "=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+\n"
				+ "Input file formats:\n"
				+ "varScan.txt has the general tab delimetted varScan layout with patient alphaNumeric (PA[0-9]{n} type): chrom	position	ref	var	normal_reads1	"
				+ "\nnormal_reads2	normal_var_freq	normal_gt	tumor_reads1	tumor_reads2	tumor_var_freq	tumor_gt	"
				+ "\nsomatic_status	variant_p_value	somatic_p_value	tumor_reads1_plus	tumor_reads1_minus	tumor_reads2_plus	"
				+ "\ntumor_reads2_minus	normal_reads1_plus	normal_reads1_minus	normal_reads2_plus	normal_reads2_minus	patient\n\n"+

				"patient_phenotype.csv is comma separated with the general layout of: Patientfactor,PHENOTYPE,Genderfactor,Locationfactor,"
				+ "\nAge,patientID,Gender,Location,BRAF.Alteration,Gene.Fusion\n\n"+
				"-hd to include column names in output; -nd for no header in output files\n\n"+
				"threshold specifies snp genotype and should be numeric (could be float)\n"
				+ "numThreads indicates the number of threads requested\n\n"
				+ "Output: varScanFileName.ped and varScanFileName.map files in the same directory. Existing files will be overwritten...\n"
				+ "\nV2.7 Developed by Fadi Hariri; fadi.hariri@mail.mcgill.ca\n");
		System.exit(0);
	}

	public static void main(String[] args) throws InterruptedException {

		/**
		File file = new File ("gwas.test.small.txt");
		File pheno = new File ("gwas.pheno.csv");
		execute (file,pheno,new File ("gwas.ped"),new File ("gwas.map"), 90, 3);
		 */
		long start=System.currentTimeMillis();

		String varscan="",pheno=""; double threshold; int numThreads;
		try{
			if (args[0].equals("--help") || args[0].equals("-h") )help();
			else{
				varscan=args[0];pheno=args[1];
				if (args[2].equals("-hd"))header=true;
				threshold=Double.parseDouble(args[3]);numThreads=Integer.parseInt(args[4]);
				String ped=varscan.split(".txt")[0]+".ped";
				String map=varscan.split(".txt")[0]+".map";
				if (threshold < 0 || numThreads <= 0){System.err.println("Numeric errors in number of threads or threshold...");System.exit(0);}
				execute (new File (varscan),new File(pheno),new File(ped),new File(map),threshold,numThreads);
			}
		}catch (ArrayIndexOutOfBoundsException e){System.err.println("missing arguments...\n --help or -h for more info.");System.exit(0);}

		long end = System.currentTimeMillis();
		System.out.println("Elapsed time: "+((double)(end-start)/60000)+" minutes");


	}

	/**
	 * Read file and parse data into proper ADTs
	 * Note File format is from varScan with modified by appending patientID as last column
	 * Will execute as a single thread
	 * @throws InterruptedException 
	 */
	private static void parser (final File file,final File pheno, final double threshold) throws InterruptedException{
		Thread parse = new Thread (new Runnable (){
			@SuppressWarnings({ "rawtypes", "unchecked" })
			public void run (){
				Scanner kb=null;
				try {
					kb = new Scanner (new BufferedReader (new FileReader (file)));
					int lines=0;//num of lines in file
					String str;
					while (kb.hasNextLine()){
						str= kb.nextLine();
						String [] att = str.split("\t");
						if (att[0].contains("chrom"))continue;//header columns
						//Attributes
						String snpName = att[0]+"-"+att[1];
						if (unwanted.contains(snpName))continue;

						String ref = att[2]+" "+att[2];
						String var= att[3]+" "+att[3];
						String mixed=att[2]+" "+att[3];
						double tumorFreq=Double.parseDouble(att[10].split("%")[0]);
						int patID = Integer.parseInt(att[23].substring(2));

						//if bad genotype syntax, then don't include
						if (ref.contains("/")||var.contains("/")||mixed.contains("/")){unwanted.add(snpName);continue;}

						//fill up data 
						//snp ref
						if (!snp_refs.containsKey(snpName))snp_refs.put(snpName, ref);
						//add patient to patient set
						all_pt.add(patID);
						//patient snps
						if (!patient_snps.containsKey(snpName)){
							//	@SuppressWarnings({ "unchecked", "rawtypes" })
							HashMap <Integer, String> arr= new HashMap ();
							if (tumorFreq >= threshold)arr.put(patID,var);
							else arr.put(patID,mixed);
							patient_snps.put(snpName,arr);
						}
						else {
							HashMap <Integer, String> arr= patient_snps.get(snpName);
							if (tumorFreq >= threshold)arr.put(patID,var);
							else arr.put(patID,mixed);
							patient_snps.remove(snpName);
							patient_snps.put(snpName,arr);
						}

						//metadata
						if (!metadata.containsKey(patID)){
							@SuppressWarnings({ })
							ArrayList <String> arr= new ArrayList ();
							arr.add(""+patID);arr.add(""+patID);arr.add("0");arr.add("0");//famID, patientID, patID, matID
							metadata.put(patID, arr);
						}
						lines++;
					}
					total=lines;

				}catch (FileNotFoundException e){System.err.println ("Error File "+file.toString()+" not found...");}
				finally {if (kb!=null)kb.close();}
			}
		});

		parse.start();
		Thread parse2 = new Thread (new Runnable (){
			public void run (){
				Scanner kb=null;
				try {
					kb = new Scanner (new BufferedReader(new FileReader (pheno)));
					int count=0;
					while (kb.hasNextLine ()){
						String str = kb.nextLine();
						count++;
						if (count == 1)continue;
						String [] att = str.split(",");
						int pt = Integer.parseInt(att[5].substring(2));
						String sex=att[6];
						String fusionFactor=att[1];

						if (!md.containsKey(pt)){
							@SuppressWarnings({ "unchecked", "rawtypes" })
							ArrayList <String> arr = new ArrayList ();
							arr.add(sex);arr.add(fusionFactor);
							md.put(pt, arr);
						}
						else {
							ArrayList <String> arr = md.get(pt);
							arr.add(sex);arr.add(fusionFactor);
							md.remove(pt);
							md.put(pt, arr);
						}
					}
				}catch (FileNotFoundException e){System.err.println ("Phenotype file not found...");}
				finally {if (kb != null)kb.close();}
			}
		});
		parse2.start();
		parse.join();
		parse2.join();

		//finilize metadata
		for (Integer patient : metadata.keySet()){metadata.get(patient).addAll(md.get(patient));}

	}


	/**
	 * Create ped and map files
	 * @throws InterruptedException 
	 */
	public static void execute (File file,File pheno, File ped,File map, double threshold, int threads) throws InterruptedException{
		System.out.println ("Parsing files...Building data structures...");
		parser (file,pheno, threshold);
		System.out.println ("Parsing complete...generating ped and map files...");
		ProgressTracker.progress=patient_snps.size();
		//create threads and pass subset of data to it
		int totalSize = patient_snps.size();
		int fractionalSize = totalSize/threads;
		int startPos =0;
		MyThread [] threading = new MyThread [(threads==1)?threads:threads+1];
		for (int i=1; i <=threads; i++){
			threading[i-1]= new MyThread ();
			setADTs (threading [i-1], startPos, startPos+fractionalSize);
			startPos+= fractionalSize;
		}

		//last thread
		if (threads > 1){
			int remainder = totalSize -fractionalSize*threads;
			threading [threads]= new MyThread ();
			setADTs (threading [threads], startPos, startPos+remainder);
		}


		for (MyThread th : threading)th.start();
		for (MyThread th : threading)th.join();


		System.out.print ("\rProgress==================================================100%\n");
		System.out.println("ped/map files completed...Writing data...");
		//done! write data

		FileWriter wr=null;
		FileWriter wr2=null;
		try {
			wr= new FileWriter (ped);
			wr2= new FileWriter (map);
			if (header){
				wr.write("fam.id\tind.id\tpat.id\tmat.id\tsex\tpheno\t");
				wr2.write("chr\tsnp\tgene.dist\tbp.pos\n");
			}
			for (String snp : all_snps){
				if (header)wr.write(snp+"\t");
				wr2.write(snp.split("-")[0]+"\t"+snp+"\t"+"0\t"+snp.split("-")[1]+"\n");
			}
			if (header)wr.write("\n");
			for (Integer pt:data.keySet()){
				for (String str : metadata.get(pt)){wr.write(str+"\t");}
				for (String str :data.get(pt)){wr.write(str+"\t");}
				wr.write("\n");
			}

		}catch (IOException e){System.err.println("File not found...");}
		finally {try {if (wr != null)wr.close();if (wr2 != null)wr2.close();}catch (IOException e){System.err.println("IO Fault...");}}
		System.out.println ("Output files: "+ped+"\t"+map);

	}


	/**
	 * copy ADTs for MyThread objects
	 */
	private static void setADTs (MyThread th, int start, int end){
		int count=0;
		for (String snp : patient_snps.keySet()){
			if (count >= start && count < end){
				//th.psnp.put(pt, patient_snp_list.get(pt));
				th.psnp_varients.put(snp, patient_snps.get(snp));
			}
			count++;
		}
		//System.out.println ("Size: "+th.psnp_varients.size());
	}

}

/**
 * Thread class for ped/map generation
 */
class MyThread extends Thread {
	/**
	@SuppressWarnings({ "unchecked", "rawtypes" })
	protected HashMap <Integer, ArrayList <String>> psnp= new HashMap ();//patient --> snp list
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	protected HashMap <String, HashMap <Integer,String>> psnp_varients= new HashMap ();//patient --> {snp-->varient}
	protected static final Lock lock = new ReentrantLock();
	//public MyThread (HashMap <Integer, ArrayList<String>> psnp,HashMap <Integer, HashMap <String,String>> psnp_varients ){this.psnp=psnp;this.psnp_varients=psnp_varients;}
	/**
	 * Currently runs in O(nm); n is snp size and m is patient size
	 */
	public void execute (){
		lock.lock();
		for (String snp : psnp_varients.keySet()){
			GWASFormatter2.all_snps.add(snp);
			for (Integer pt : GWASFormatter2.all_pt){
				if (psnp_varients.get(snp).keySet().contains(pt)){//snp exists in patient data then add proper snp else add ref genotype
					if (!GWASFormatter2.data.containsKey(pt)){
						@SuppressWarnings({ "unchecked", "rawtypes" })
						ArrayList <String> arr= new ArrayList ();
						arr.add(psnp_varients.get(snp).get(pt));//get that patient genotype for the specific snp
						GWASFormatter2.data.put(pt, arr);

					}
					else {
						ArrayList <String> arr= GWASFormatter2.data.get(pt);
						arr.add(psnp_varients.get(snp).get(pt));
						GWASFormatter2.data.remove(pt);
						GWASFormatter2.data.put(pt, arr);
					}
				}
				else {//add ref genotype
					if (!GWASFormatter2.data.containsKey(pt)){
						@SuppressWarnings({ "unchecked", "rawtypes" })
						ArrayList <String> arr= new ArrayList ();
						arr.add(GWASFormatter2.snp_refs.get(snp));
						GWASFormatter2.data.put(pt, arr);
					}
					else {
						ArrayList <String> arr= GWASFormatter2.data.get(pt);
						arr.add(GWASFormatter2.snp_refs.get(snp));
						GWASFormatter2.data.remove(pt);
						GWASFormatter2.data.put(pt, arr);
					}
				}
			}
			ProgressTracker.progressBar();
		}
		//System.out.println ("Data size:"+GWASFormatter2.data.size());
		lock.unlock();
	}

	public void run (){
		this.execute();
	}
}
/**
 * Progress Tracker
 */
class ProgressTracker {
	static int progress =GWASFormatter2.size;
	static int totalCount = progress/100;
	public static void progressBar (){
		if (((double)progress/totalCount)%2 == 0){
			//progString += "=";
			String out="\rProgress";
			int eq=(100-progress/totalCount)/2;//= for each 2%
			for (int i=1; i <= eq;i++)out+= "=";
			int numSpaces=50-eq;
			for (int i=1; i <= numSpaces;i++)out+= " ";
			System.out.print (out+(eq*2)+" %");
		}
		progress--;
	}
}