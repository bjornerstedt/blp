tab = array2table([12,23;35,3]);
tab.Properties.RowNames = {'A','B'}
tab2 = array2table([12,23;35,3;35,3]);
tab2.Properties.RowNames = {'A','B','D'}
et = ExcelTable('test.xlsx');
et.setSheet('huto');
et.write(tab,  'heading', 'Testing', 'rowheading', 'test');
et.write(tab2,  'below', false);
et.write(tab,  'rowheading', 'test', 'heading', 'Testing');

et.setSheet('josefin');
et.write(tab2, 'below', false);
et.write(tab2);
et.setSheet('test');

et.write(tab, 'sheet', 'test', 'below', false);
et.write(tab2, 'sheet', 'test');
et.write(tab2, 'sheet', 'test', 'below', false);
