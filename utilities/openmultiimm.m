function out = openmultiimm(varargin)
% ---
% --- OPENMULTIIMM Open *.imm multifiles
% ---
% --- USAGE : OUT = OPENMULTIIMM(FILENAME,INDEX)
% ---
% --- Input Argument:
% --- FILENAME       : image file name
% --- INDEX          : image index number
% --- IMAGESTARTBYTE : previously by INDEXCOMPRESSEDMULTIIMM determined
% ---                  starting byte positition of image #index
% ---
% --- Output Argument:
% --- OUT : structure containing header & image. Usually image is an
% ---       intensity image with 0 for black, saturation point for white
% ---       (65535 for images saved by Yorick).
% ---
% --- Zhang Jiang & Michael Sprung
% --- $Revision: 1.0 $Date: 2004/12/09 $
% --- $Revision: 1.1 $Date: 2005/12/01 $ change to use updated imm header
% --- $Revision: 1.2 $Date: 2006/03/06 $ change to use updated imm header
% --- $Revision: 2.0 $Date: 2006/09/06 $ change to open either compressed
% ---                                    or uncompresses *.imm multifiles
% ---
% --- $Revision: 2.1 $Date: 2015/08/03 : cleaned up quite a bit, removed
% undefined headers, changed the file index so for a given imm file, the
% first frame# is always 1 and not depending on the file name suffix
% =========================================================================
% --- input file preparation
% =========================================================================
file     = varargin{1}                                                     ;
immIndex = varargin{2}                                                     ;

if ( nargin == 3 && ~isempty(varargin{3}) == 1 )
    imageStart = varargin{3}                                           ; % image starting position is already provided
else
    imageStart = [];
end
% =========================================================================
% --- open file
% =========================================================================
[fid,message] = fopen(file)                                                ;        
if ( fid == -1 )                                                             % return if open fails
    uiwait(msgbox(message,'File Open Error','error','modal'))              ;
    return                                                                 ;
end
% =========================================================================
% --- check for mode, compression & immversion
% =========================================================================
fseek(fid,0,'bof')                                                         ;
modeflag        = fread(fid,1 ,'int')                                      ;
fseek(fid,4,'bof')                                                         ;
compressionflag = fread(fid,1 ,'int')                                      ;
fseek(fid,616,'bof')                                                       ;
immversionflag  = fread(fid,1 ,'int')                                      ;
% --- check for multifile format
if ( immversionflag >= 11 && modeflag ~= 2 )
    uiwait(msgbox(message,'Not a Multifile','error','modal'))              ;    
    fclose(fid)                                                            ;
    return                                                                 ;
end
% --- check for compression
compression = 0                                                            ; % assume uncompressed files
if ( compressionflag == 6 && immversionflag >= 11 )                      % condition for compressed data from Imageserver program
    compression = 1                                                        ;
end

% =========================================================================
% --- determine the start position of the wanted image
% =========================================================================
fseek(fid,108,'bof')                                                       ;
rows = fread(fid,1 ,'int')                                                 ;
fseek(fid,112,'bof')                                                       ;
cols = fread(fid,1 ,'int')                                                 ;
fseek(fid,116,'bof')                                                       ;
bytes = fread(fid,1 ,'int')                                                ;

if ( (compression == 0 ) && isempty(imageStart) )
    if ( immversionflag >= 11 )
        immSize = bytes * rows * cols + 1024                               ; % 'bytes' per pixel
    else
        immSize = 2 * rows * cols + 1024                                   ; % always 2 bytes per pixel in older versions
    end
    imageStart = (immIndex-1)*immSize                                      ;
    
elseif ( (compression == 1 ) && isempty(imageStart) )
    imageStart      = 0                                                ;
    for k = 1 : immIndex - 1
        fseek(fid,imageStart+152,'bof')                                ;
        dlen = fread(fid,1 ,'int')                                     ; % 'dlen' in this case contains the # of exposed pixels
        imageStart = imageStart + 1024 + dlen * ( 4 + bytes )          ;
    end
end

% =========================================================================
% --- load header of given image index
% =========================================================================
fseek(fid,imageStart,'bof')                                                ;
header = cell(53,2)                                                        ; % initialize '.header' structure
header(1,:)  = {'mode',          fread(fid,1 ,'int')}                      ; % starting byte   0 : length  4 bytes
header(2,:)  = {'compression',   fread(fid,1 ,'int')}                      ; % starting byte   4 : length  4 bytes
header(3,:)  = {'date',          fread(fid,32,'uchar')'}                   ; % starting byte   8 : length 32 bytes
% header(4,:)  = {'prefix',        fread(fid,16,'uchar')'}                   ; % starting byte  40 : length 16 bytes
% header(5,:)  = {'number',        fread(fid,1 ,'int')}                      ; % starting byte  56 : length  4 bytes
% header(6,:)  = {'suffix',        fread(fid,16,'uchar')'}                   ; % starting byte  60 : length 16 bytes
% header(7,:)  = {'monitor',       fread(fid,1 ,'int')}                      ; % starting byte  76 : length  4 bytes
% header(8,:)  = {'shutter',       fread(fid,1 ,'int')}                      ; % starting byte  80 : length  4 bytes
fseek(fid,imageStart+84,'bof');
header(9,:)  = {'row_beg',       fread(fid,1 ,'int')}                      ; % starting byte  84 : length  4 bytes
header(10,:) = {'row_end',       fread(fid,1 ,'int')-1}                    ; % starting byte  88 : length  4 bytes
header(11,:) = {'col_beg',       fread(fid,1 ,'int')}                      ; % starting byte  92 : length  4 bytes
header(12,:) = {'col_end',       fread(fid,1 ,'int')-1}                    ; % starting byte  96 : length  4 bytes
header(13,:) = {'row_bin',       fread(fid,1 ,'int')}                      ; % starting byte 100 : length  4 bytes
header(14,:) = {'col_bin',       fread(fid,1 ,'int')}                      ; % starting byte 104 : length  4 bytes
header(15,:) = {'rows',          fread(fid,1 ,'int')}                      ; % starting byte 108 : length  4 bytes
header(16,:) = {'cols',          fread(fid,1 ,'int')}                      ; % starting byte 112 : length  4 bytes
header(17,:) = {'bytes',         fread(fid,1 ,'int')}                      ; % starting byte 116 : length  4 bytes
header(18,:) = {'kinetics',      fread(fid,1 ,'int')}                      ; % starting byte 120 : length  4 bytes
header(19,:) = {'kinwinsize',    fread(fid,1 ,'int')}                      ; % starting byte 124 : length  4 bytes
header(20,:) = {'elapsed',       fread(fid,1 ,'double')}                   ; % starting byte 128 : length  8 bytes
header(21,:) = {'preset',        fread(fid,1 ,'double')}                   ; % starting byte 136 : length  8 bytes
% header(22,:) = {'topup',         fread(fid,1 ,'int')}                      ; % starting byte 144 : length  4 bytes
% header(23,:) = {'inject',        fread(fid,1 ,'int')}                      ; % starting byte 148 : length  4 bytes
fseek(fid,imageStart+152,'bof');
header(24,:) = {'dlen',          fread(fid,1 ,'int')}                      ; % starting byte 152 : length  4 bytes
header(25,:) = {'roi_number',    fread(fid,1 ,'int')}                      ; % starting byte 156 : length  4 bytes
header(26,:) = {'buffer_number', fread(fid,1 ,'uint32')}                   ; % starting byte 160 : length  4 bytes
header(27,:) = {'systick',       fread(fid,1 ,'uint32')}                   ; % starting byte 164 : length  4 bytes
% --- shifted header positions as of 20060306
% header(28,:) = {'pv1',           fread(fid,40,'uchar')}                    ; % starting byte 168 : length 40 bytes
% header(29,:) = {'pv1VAL',        fread(fid,1 ,'float')}                    ; % starting byte 208 : length  4 bytes
% header(30,:) = {'pv2',           fread(fid,40,'uchar')}                    ; % starting byte 212 : length 40 bytes
% header(31,:) = {'pv2VAL',        fread(fid,1 ,'float')}                    ; % starting byte 252 : length  4 bytes
% header(32,:) = {'pv3',           fread(fid,40,'uchar')}                    ; % starting byte 256 : length 40 bytes
% header(33,:) = {'pv3VAL',        fread(fid,1 ,'float')}                    ; % starting byte 296 : length  4 bytes
% header(34,:) = {'pv4',           fread(fid,40,'uchar')}                    ; % starting byte 300 : length 40 bytes
% header(35,:) = {'pv4VAL',        fread(fid,1 ,'float')}                    ; % starting byte 340 : length  4 bytes
% header(36,:) = {'pv5',           fread(fid,40,'uchar')}                    ; % starting byte 344 : length 40 bytes
% header(37,:) = {'pv5VAL',        fread(fid,1 ,'float')}                    ; % starting byte 384 : length  4 bytes
% header(38,:) = {'pv6',           fread(fid,40,'uchar')}                    ; % starting byte 388 : length 40 bytes
% header(39,:) = {'pv6VAL',        fread(fid,1 ,'float')}                    ; % starting byte 428 : length  4 bytes
% header(40,:) = {'pv7',           fread(fid,40,'uchar')}                    ; % starting byte 432 : length 40 bytes
% header(41,:) = {'pv7VAL',        fread(fid,1 ,'float')}                    ; % starting byte 472 : length  4 bytes
% header(42,:) = {'pv8',           fread(fid,40,'uchar')}                    ; % starting byte 476 : length 40 bytes
% header(43,:) = {'pv8VAL',        fread(fid,1 ,'float')}                    ; % starting byte 516 : length  4 bytes
% header(44,:) = {'pv9',           fread(fid,40,'uchar')}                    ; % starting byte 520 : length 40 bytes
% header(45,:) = {'pv9VAL',        fread(fid,1 ,'float')}                    ; % starting byte 560 : length  4 bytes
% header(46,:) = {'pv10',          fread(fid,40,'uchar')}                    ; % starting byte 564 : length 40 bytes
% header(47,:) = {'pv10VAL',       fread(fid,1 ,'float')}                    ; % starting byte 604 : length  4 bytes
% --- new fields (added 10/2006)
% header(48,:) = {'imageserver',   fread(fid,1 ,'float')}                    ; % starting byte 608 : length  4 bytes
% header(49,:) = {'CPUspeed',      fread(fid,1 ,'float')}                    ; % starting byte 612 : length  4 bytes
fseek(fid,imageStart+616,'bof');
header(50,:) = {'immversion',    fread(fid,1 ,'int')}                      ; % starting byte 616 : length  4 bytes
header(51,:) = {'corecotick',    fread(fid,1 ,'uint32')}                   ; % starting byte 620 : length  4 bytes
% header(52,:) = {'cameratype',    fread(fid,1 ,'int')}                      ; % starting byte 624 : length  4 bytes
% header(53,:) = {'threshhold',    fread(fid,1 ,'float')}                    ; % starting byte 628 : length  4 bytes


% =========================================================================
% --- store header to output structure
% =========================================================================
out.header = header                                                        ;


% =========================================================================
% --- store image to output structure
% =========================================================================
if ( compression == 0 )

    fseek(fid,imageStart+1024,'bof')                                       ;
    if ( immversionflag >= 11 )
        if ( bytes == 2 )
            out.imm = fread(fid,[header{16,2},header{15,2}],'uint16=>single')' ;
        elseif ( bytes == 4 )
            out.imm = fread(fid,[header{16,2},header{15,2}],'uint32=>single')' ;
        end
    else
        out.imm = fread(fid,[header{16,2},header{15,2}],'uint16=>single')' ;
    end

elseif ( compression == 1 )
    
    if ( compressionflag == 6 && immversionflag >= 11 )
        PixelNumber = floor(double(header{24,2}))                          ; % for compressed mode 'dlen' is the # of exposed pixels
        fseek(fid ,imageStart + 1024,'bof')                                ; % go to the next data start of the compressed file
        PixelIndex = fread(fid,PixelNumber,'uint32') + 1                   ; % read the full index block (correction to start at 1)
        if ( bytes == 2 )
            PixelValue = fread(fid,PixelNumber,'uint16=>single')           ; % read the full value block
        elseif ( bytes == 4 )
            PixelValue = fread(fid,PixelNumber,'uint32=>single')           ; % read the full value block
        end
        % --- create uncompressed data
        UCData = zeros(header{16,2},header{15,2},'single')                 ;
        UCData(PixelIndex) = PixelValue                                    ;
        % --- store image data
        out.imm = UCData'                                                  ;
    end
    
    
end


% =========================================================================
% --- close fid and return;
% =========================================================================
fclose(fid)                                                                ;


% ---
% EOF
