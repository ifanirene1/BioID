//This macro will separate a 3 channel image that is stained
dir = getDirectory("Select Source Directory");
//define the output directory and create it 
saveDir = dir+"Output"; 
File.makeDirectory(saveDir);
fs = File.separator;
//enter batch mode to speed up processing
//setBatchMode(true);
//generate file list for batch processing 
list = getFileList(dir);
Array.sort(list);
//loop for all images in the selected directory to count and 
//calculate Fluorecent intensites
for(i=0;i<list.length;i++)
{ filename = dir +list[i];
if (endsWith(filename, "czi")) 
{
		run("Bio-Formats", "open=["+filename+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_1");
		nameStore = getTitle();
		selectWindow(nameStore);
		makeRectangle(0, 190, 1020, 420);
		run("Crop");
           run("Subtract Background...", "rolling=10 light disable stack");
		run("Enhance Contrast...", "saturated=1 equalize process_all use");
		run("Find Edges", "stack");
		//run("Sharpen", "stack");
		run("Gaussian Blur...", "sigma=2 stack");
		setAutoThreshold("Li dark stack");
		run("Set Measurements...", "area area_fraction limit display redirect=None decimal=3");
		run("Analyze Particles...", "size=0-Infinity show=Masks clear stack");
		Mask = "Mask of "+nameStore;
		selectWindow(Mask);
		run("Invert", "stack");
		run("Analyze Particles...", "size=5000-Infinity show=Nothing clear include summarize stack");
		IJ.renameResults("Results");
		wait(1000);
		saveAs("Results", saveDir+fs+nameStore+".csv");
		//close all open images
		run("Close All");
		close("Results");
}
}
waitForUser("Ding Dong!","That was fun! What's next?")