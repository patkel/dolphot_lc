************
Analyze
************


Refer to the /example/analyze_data.py file in order to garner results. Add the following code to this file.

Declare the type of object to be analyzed and pixel coordinates of the object. ::

    objs = [['SN', 567.9625526342381, 176.49369319498047]]

Finally, run: ::

    imroot = os.getcwd()                                # example folder / base directory
    imdir_prepped = f'{imroot}/diffs/'
    a = dolphot_output( imdir_prepped + 'singlestar')
    for objname,x,y in objs:
        a.info_single_object(x, y, objname, 'orig', settings=False, verbose=False)
