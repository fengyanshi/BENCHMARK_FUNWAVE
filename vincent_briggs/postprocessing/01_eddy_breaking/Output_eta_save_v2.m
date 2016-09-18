clear all; clc; close all;

inputpath='../result_01/';
fn=load('./fname.in');

for k=1:1:length( fn )
    
        fname=[ inputpath, 'eta_', num2str( fn(k), '%05d' ) ];
        eta = load( fname );

        outname=[ inputpath, 'eta_', num2str( fn(k), '%05d' ), '.mat' ];
        save( outname, 'eta' );
        clear eta;
        
%         delete( fname )
        clear fname outname
 
end