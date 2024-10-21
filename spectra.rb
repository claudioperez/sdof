
class Spectra < ApplicationRecord
  has_many :searches
  has_many :unscaled_searches

  DIRECTORY = File.join("public", "user_defined_spectrum")

  HUMANIZED_COLLUMNS = {
      :NGAInputData_M => "Magnitude",
      :NGAInputData_Vs30 => "Vs30",
      :NGAInputData_RRUP => "RRUP",
      :NGAInputData_DIP => "DIP",
      :NGAInputData_ZTOR => "ZTOR",
      :NGAInputData_W => "Width",
      :NGAInputData_RJB => "RJB",
      :NGAInputData_Z25 => "Z.25",
      :NGAInputData_Z10 => "Z.10",
      :NGAInputData_Rx => "Rx",
      :NGAInputData_HWflag => "Hanging Wall",
      :NGAMeanFlag => "GMPE average",
      :NGAInputData_CMFlag => "CMS Flag",
      :NGAInputData_Vs30_class => "Vs30 (estimated)"
  }

 # GV Comment Start
   # AS BA CB CY ID

   validates :NGAInputData_M, presence: true, :if => Proc.new { |spectra| spectra.NGAInputData_NGAModelSelection == 1}
   validates :menu_Mechanism, presence: true, :if => Proc.new { |spectra| spectra.NGAInputData_NGAModelSelection == 1}
   validates :NGAInputData_Vs30, presence: true, :if => Proc.new { |spectra| spectra.NGAInputData_NGAModelSelection == 1}
   validates :NGAInputData_HWflag, presence: true, :if => Proc.new { |spectra| spectra.NGAInputData_NGAModelSelection == 1}

   validates :NGAMeanFlag,presence: true,  :if => Proc.new { |spectra| spectra.NGAInputData_NGAModelSelection == 1}
   validates :NGAInputData_CMFlag, presence: true, :if => Proc.new { |spectra| spectra.NGAInputData_NGAModelSelection == 1}


   validates :NGAInputData_Vs30,presence: true,  :if => Proc.new { |spectra| spectra.validate_fields_vs30?(0,0,0,0,1, spectra.NGAInputData_checkboxNGA, spectra.NGAInputData_Vs30) }


   validates :NGAInputData_Tdesign, presence: true, :if => Proc.new { |spectra| spectra.validate_fields_Tdesign?(1,1,1,1,1, spectra.NGAInputData_CMFlag, spectra.NGAInputData_Tdesign) }
   validates :NGAInputData_EpsTdesign,presence: true,  :if => Proc.new { |spectra| spectra.validate_fields_EpsTdesign?(1,1,1,1,1, spectra.NGAInputData_CMFlag) }
   # AS    CB CY ID  - corresponds to the as=1,ba=0,cb=1,cy=1,id=1 passed to validate
   validates :NGAInputData_RRUP,presence: true,  :if => Proc.new { |spectra| spectra.validate_fields?(1,0,1,1,1, spectra.NGAInputData_checkboxNGA) }
  validates :NGAInputData_DIP, presence: true,  :if => Proc.new { |spectra| spectra.validate_fields?(1,0,1,1,0, spectra.NGAInputData_checkboxNGA) }
  validates :NGAInputData_ZTOR,presence: true,  :if => Proc.new { |spectra| spectra.validate_fields?(1,0,1,1,0, spectra.NGAInputData_checkboxNGA) }
  validates :NGAInputData_W, presence: true,    :if => Proc.new { |spectra| spectra.validate_fields?(1,0,1,0,0, spectra.NGAInputData_checkboxNGA) }
  validates :NGAInputData_RJB,presence: true,   :if => Proc.new { |spectra| spectra.validate_fields?(1,1,1,1,0, spectra.NGAInputData_checkboxNGA) }
  validates :NGAInputData_Region, presence: true,  :if => Proc.new { |spectra| spectra.validate_fields?(1,1,1,1,0, spectra.NGAInputData_checkboxNGA) }
  validates :NGAInputData_Z25,presence: true,   :if => Proc.new { |spectra| spectra.validate_fields?(0,0,1,0,0, spectra.NGAInputData_checkboxNGA) }
  validates :NGAInputData_Z10, presence: true,  :if => Proc.new { |spectra| spectra.validate_fields?(1,1,0,1,0, spectra.NGAInputData_checkboxNGA) }
  validates :NGAInputData_ZHYP, presence: true,  :if => Proc.new { |spectra| spectra.validate_fields?(0,0,1,0,0, spectra.NGAInputData_checkboxNGA) }
  validates :NGAInputData_Rx,presence: true,    :if => Proc.new { |spectra| spectra.validate_fields?(1,0,1,1,0, spectra.NGAInputData_checkboxNGA) }
  validates :NGAInputData_RY0,presence: true,    :if => Proc.new { |spectra| spectra.validate_fields?(1,0,0,0,0, spectra.NGAInputData_checkboxNGA) }
  validates :NGAInputData_Vs30_class,presence: true,  :if => Proc.new { |spectra| spectra.validate_fields?(1,0,0,1,0, spectra.NGAInputData_checkboxNGA) }
  validates :w_ask14,presence: true,  :if => Proc.new { |spectra| spectra.validate_fields?(1,0,0,0,0, spectra.NGAInputData_checkboxNGA) }
  validates :w_bssa14, presence: true, :if => Proc.new { |spectra| spectra.validate_fields?(0,1,0,0,0, spectra.NGAInputData_checkboxNGA) }
  validates :w_cb14,presence: true,  :if => Proc.new { |spectra| spectra.validate_fields?(0,0,1,0,0, spectra.NGAInputData_checkboxNGA) }
  validates :w_cy14,presence: true,  :if => Proc.new { |spectra| spectra.validate_fields?(0,0,0,1,0, spectra.NGAInputData_checkboxNGA) }
  validates :w_i14, presence: true, :if => Proc.new { |spectra| spectra.validate_fields?(0,0,0,0,1, spectra.NGAInputData_checkboxNGA) }
 # GV Comment End


 # GV Comment Start
 validates :filename, presence: true, :if => Proc.new { |spectra| spectra.NGAInputData_NGAModelSelection == 0}
 # GV Comment End


  def CreateSpectraImage(params,targetSpectrumArray)
    puts "***********************"
    puts "***********************"
    puts " "
    puts 'Start in spectras model -> CreateSpectraImage'
    puts " "
    puts "***********************"
    puts "***********************"

    puts "targetSpectrumArray = "+targetSpectrumArray.to_s

    # startTime=Time.now
    thisSpectraID=self.id

    #@spectras = Spectra.all

    nGAModelSelection = self.NGAInputData_NGAModelSelection


    # prepared_stmt ="delete _spectra_selectedtargetspectrum.* from _spectra_selectedtargetspectrum where (((_spectra_selectedtargetspectrum.spectraid) =#{thisSpectraID})); "
    #    ActiveRecord::Base.connection.execute(prepared_stmt)
    #      ActiveRecord::Base.connection_handler.clear_active_connections!



    targetPeriodList = targetSpectrumArray["targetPeriodList"]
    targetpSaList = targetSpectrumArray["targetpSaList"]
    targetpSaArray = targetSpectrumArray["targetpSaArray"]
    allGMPEplotList = targetSpectrumArray["allGMPEplotList"]
    allGMPEplotDataList = targetSpectrumArray["allGMPEplotDataList"]


    periodStartCounter = 0
    targetSpectrumPlotList = ""
    targetPeriodList.each {|thisT|
      unless thisT>0 then next end
      thispSa = targetpSaArray[thisT]
      #    prepared_stmt = "  INSERT INTO _spectra_selectedtargetspectrum (spectraid,tperiod,psatarget,tsa)  VALUES (#{thisSpectraID},#{thisT},#{thispSa},\"#{thisT}\,#{thispSa}\")"
      #   ActiveRecord::Base.connection.execute(prepared_stmt)
      # ActiveRecord::Base.connection_handler.clear_active_connections!
      if periodStartCounter == 0
        targetSpectrumPlotList = "[[#{thisT},#{thispSa}]"
        periodStartCounter = 1
      else
        targetSpectrumPlotList = "#{targetSpectrumPlotList},[#{thisT},#{thispSa}]"
      end
    }
    targetSpectrumPlotList = "#{targetSpectrumPlotList}]"
    targetSpectrumPlotDataList = "{linePattern: 'solid', lineWidth: 3, color: 'red', label: 'Target Spectrum'}"

    if (nGAModelSelection == 1 )

      targetSpectrumPlotList = "#{allGMPEplotList},#{targetSpectrumPlotList}"
      targetSpectrumPlotDataList = "#{allGMPEplotDataList},#{targetSpectrumPlotDataList}"


    end



    unique_token = getRandomString
    self.unique_token = unique_token


    commandList = [
        "update spectras set spectras.targetSpectrumData = null where spectras.id=#{thisSpectraID}; ",
        "update spectras set spectras.TargetSpectrumSeriesProps = null where spectras.id=#{thisSpectraID}; ",
        "update spectras set targetSpectrumData = \"#{targetSpectrumPlotList}\" where spectras.id=#{thisSpectraID}; ",
        "update spectras set TargetSpectrumSeriesProps = \"#{targetSpectrumPlotDataList}\" where spectras.id=#{thisSpectraID}; ",
        "update spectras set unique_token = \"#{unique_token}\" where spectras.id=#{thisSpectraID}; ",
    ]

    commandList.each do |prepared_stmt|
      (ActiveRecord::Base.connection.execute(prepared_stmt))
    end
    ActiveRecord::Base.connection_handler.clear_active_connections!




    #     calcTime = Time.now-startTime
  end




  def CreateTargetSpectra_Mapped(params)

    startTime=Time.now
    thisSpectraID=self.id

    #@spectras = Spectra.all
    puts "in CreateTargetSpectra_Mapped "
    puts self.id

    nGAModelSelection = self.NGAInputData_NGAModelSelection

    # debugger
    targetSpectrumArray = {}
    targetPeriodList = []
    targetpSaList = []
    targetpSaArray = {}




    prepared_stmt = "SELECT id,created_at,updated_at,NGAInputData_NGAModelSelection,user_id,ngainputdata_sds,ngainputdata_sd1,ngainputdata_tl
                        FROM spectras where (id = #{thisSpectraID});"
    result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
    ActiveRecord::Base.connection_handler.clear_active_connections!

    id,created_at,updated_at,ngamodelselection,user_id,sds,sd1,tL = result_set.first

    t0 = 0.2*sd1/sds
    t1 = sd1/sds


    prepared_stmt = "SELECT tperiod FROM z_recordspectraperiods order by tperiod"
    result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
    ActiveRecord::Base.connection_handler.clear_active_connections!

    result_set.each {|thisRow|
      thisT = thisRow[0]

      if thisT < t0
        thispSa = sds*(0.4+0.6*thisT/t0)
      elsif thisT <=t1
        thispSa = sds
      elsif thisT <= tL
        thispSa = sd1/thisT
      else
        thispSa = sd1*tL/thisT**2
      end

      targetPeriodList << thisT
      targetpSaList << thispSa
      targetpSaArray[thisT] = thispSa
    }


    targetSpectrumArray["targetPeriodList"] = targetPeriodList
    targetSpectrumArray["targetpSaList"] = targetpSaList
    targetSpectrumArray["targetpSaArray"] = targetpSaArray

    return  targetSpectrumArray
  end




  def CreateTargetSpectra_UserDefined(params)

    startTime=Time.now
    thisSpectraID=self.id

    #@spectras = Spectra.all
    logger.info "in CreateTargetSpectra_UserDefined "
    # logger.info self.id

    nGAModelSelection = self.NGAInputData_NGAModelSelection

    # debugger

    directory = "public/upload"
    filepath = directory+"/"+self.filename.to_s
    logger.info 'filepath'
    logger.info filepath
    # file = filepath

    targetSpectrumArray = {}
    targetPeriodList = []
    targetpSaList = []
    targetpSaArray = {}
    if (nGAModelSelection == 0 )
      if (filename!= nil and filename!="")

        if File.file?(filepath) == false
          filepath = filepath.gsub("+", " ")
        end
        if File.file?(filepath) == false
          return -1
        end
        logger.info 'start'
        # 150914 silviamazzoni -- read the file this way so that you may deal with mac
        lineCounter = 0
        File.open(filepath).each(sep="\r\n") do |line|

          logger.info 'in first IF rn'

          # logger.info 'line'
          # logger.info line
          lineCounter += 1
          # next if lineCounter<4
          next if line.length == 0

          parsed_points = line.strip.split(',')
          # logger.info 'parsed_points'
          # logger.info parsed_points

          next if parsed_points.length != 2

          t = parsed_points[0].to_f
          pSa = parsed_points[1].to_f
          if (t - 0.0<0.0001)
            t=0.00001
          end
          if (t>0 and pSa>0)
            # logger.info 'both greater than zero'
            targetPeriodList << t
            targetpSaList << pSa
            targetpSaArray[t]=pSa
          end

        end

        logger.info 'pc type rn'
        logger.info lineCounter


        if lineCounter<2
          lineCounter = 0

          logger.info 'before second IF'
          File.open(filepath).each(sep="\r") do |line|
            # logger.info 'in second IF'

            # logger.info 'line'
            # logger.info line

            lineCounter += 1
            # next if lineCounter<4
            next if line.length == 0

            parsed_points = line.strip.split(',')

            # logger.info 'parsed_points'
            # logger.info parsed_points

            # logger.info 'parsed_points_length'
            # logger.info parsed_points.length

            next if parsed_points.length != 2

            t = parsed_points[0].to_f
            pSa = parsed_points[1].to_f
            # logger.info 't'
            # logger.info t
            # logger.info 'pSa'
            # logger.info pSa
            if (t - 0.0<0.0001)
              t=0.00001
            end
            if (t>0 and pSa>0)
              # logger.info 'both greater than zero'
              targetPeriodList << t
              targetpSaList << pSa
              targetpSaArray[t]=pSa
            end

          end

          logger.info 'mac Type one r'
          logger.info lineCounter
        end


        if lineCounter<2
          lineCounter = 0

          logger.info 'before third IF'
          File.open(filepath).each(sep="\n\n") do |line|
            # logger.info 'in third IF'

            # logger.info 'line'
            # logger.info line

            lineCounter += 1
            # next if lineCounter<4
            next if line.length == 0

            parsed_points = line.strip.split(',')

            # logger.info 'parsed_points'
            # logger.info parsed_points

            # logger.info 'parsed_points_length'
            # logger.info parsed_points.length

            next if parsed_points.length != 2

            t = parsed_points[0].to_f
            pSa = parsed_points[1].to_f
            # logger.info 't'
            # logger.info t
            # logger.info 'pSa'
            # logger.info pSa
            if (t - 0.0<0.0001)
              t=0.00001
            end
            if (t>0 and pSa>0)
              # logger.info 'both greater than zero'
              targetPeriodList << t
              targetpSaList << pSa
              targetpSaArray[t]=pSa
            end

          end

          logger.info 'mac Type two nn'
          logger.info lineCounter
        end

        logger.info 'done with ifs'

      end
    end

    logger.info 'targetPeriodList'
    logger.info targetPeriodList
    logger.info 'targetpSaList'
    logger.info targetpSaList
    logger.info 'targetpSaArray'
    logger.info targetpSaArray

    targetSpectrumArray["targetPeriodList"] = targetPeriodList
    targetSpectrumArray["targetpSaList"] = targetpSaList
    targetSpectrumArray["targetpSaArray"] = targetpSaArray

    logger.info 'done with userdefined target spectrum'
    return  targetSpectrumArray
  end






  def CreateTargetSpectra_GMPEs(thisSpectraID,params,models)

     puts 'MODEL = CreateTargetSpectra_GMPEs'
     puts "thisSpectraID = #{thisSpectraID}"
     puts "parmas = "+params.to_s
     puts "models = "+models.to_s


    targetSpectrumGMPEArray = {}
    gmpeList = []
    gmpeWeightSum = 0
    thisWeight = 0


    modelGMPEmapArray = {}
    modelGMPEmapArray["AS"] = "ask14"
    modelGMPEmapArray["CB"] = "cb14"
    modelGMPEmapArray["BA"] = "bssa14"
    modelGMPEmapArray["CY"] = "cy14"
    modelGMPEmapArray["ID"] = "i14"

    keyGMPEColorArray = {}
    keyGMPEColorArray["ask14"] = "blue"
    keyGMPEColorArray["cb14"] = "magenta"
    keyGMPEColorArray["bssa14"] = "black"
    keyGMPEColorArray["cy14"] = "cyan"
    keyGMPEColorArray["i14"] = "grey"

    keyGMPELabelArray = {}
    keyGMPELabelArray["ask14"] = "Abrahamson-Silva-Kamai 2014 median"
    keyGMPELabelArray["cb14"] = "Campbell-Bozorgnia 2014 median"
    keyGMPELabelArray["bssa14"] = "Boore-Stewart-Seyhan-Atkinson 2014 median"
    keyGMPELabelArray["cy14"] = "Chiou-Youngs 2014 median"
    keyGMPELabelArray["i14"] = "Idriss 2014 median"




# ----------------------------------------
    prepared_stmt = "SELECT id,created_at,updated_at,NGAInputData_NGAModelSelection,NGAInputData_checkboxNGA,NGAInputData_Nsigma,menu_Mechanism,NGAInputData_M,NGAInputData_RRUP,NGAInputData_RJB,NGAInputData_Rx,NGAInputData_HWflag,NGAInputData_ZTOR,NGAInputData_ZHYP,NGAInputData_DIP,NGAInputData_Vs30,NGAInputData_Vs30_class,NGAInputData_Z10,NGAInputData_Z25,NGAInputData_W,NGAInputData_RY0,NGAInputData_FAS,NGAInputData_Region,NGAInputData_deltaDPP,NGAInputData_CMFlag,NGAInputData_EpsTdesign,NGAInputData_Tdesign,NGAInputData_Sds,NGAInputData_Sd1,NGAInputData_TL,temp,NGAMeanFlag,varargin,unique_token,grid_on,xyScale,filename,NGAInputData_ModelCheck,NGAInputData_DampingRatio,user_id,w_ask14,w_bssa14,w_cb14,w_cy14,w_i14
                        FROM spectras where (id = #{thisSpectraID});"
    result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
    ActiveRecord::Base.connection_handler.clear_active_connections!

    id,created_at,updated_at,ngamodelselection,checkboxnga,nsigma,menu_mechanism,magnitude,rrup,rjb,rx,hwflag,ztor,zhyp,dip,vs30,vs30_class,z10,z25,w,ry0,fas,region,deltadpp,cmflag,epstdesign,tdesign,sds,sd1,tl,temp,ngameanflag,varargin,unique_token,grid_on,xyscale,filename,modelcheck,dampingratio,user_id,w_ask14,w_bssa14,w_cb14,w_cy14,w_i14  = result_set.first

    pi = 4 * Math.atan(1)
    u = 0;  # since we always define a mechanism, u is always zero, never undefined
    fhw = hwflag
    fas = 0
    vs30class=vs30_class


    if (z10 == nil)
      z10 = 999
    else
      if (z10[0..1].to_f >0)
        z10 = z10.to_f
        # z1_ask14 = z10.to_f
      else
        z10 = 999
      end
    end

    keyNGAMeanArray={}
    keyNGAMeanArray[0] = "Geometric"
    keyNGAMeanArray[1] = "Arithmetic"
    ngameantype = keyNGAMeanArray[ngameanflag]


    keyMenuMechanismArray={}
    keyMenuMechanismArray[1]="Strike Slip"
    keyMenuMechanismArray[2]="Normal/Oblique"
    keyMenuMechanismArray[3]="Reverse/Oblique"
    faulttype = keyMenuMechanismArray[menu_mechanism.to_i]


    if faulttype == "Normal/Oblique"
      frv,fnm,fss = 0,1,0
    elsif faulttype == "Reverse/Oblique"
      frv,fnm,fss = 1,0,0
    else
      frv,fnm,fss = 0,0,1
    end
    fidriss = frv


    ztord_cb14=frv*([2.704-1.226*[magnitude-5.849,0].max,0].max)**22 + (1-frv)*([2.673-1.136*[magnitude-4.97,0].max,0].max)**2
    if ztord_cb14>20
      ztord_cb14 = 0
    end


    if magnitude<6.75
      term1 = -4.317+0.984*magnitude
    else
      term1 = 2.325
    end
    if dip<40
      term2 = 0.0445*(dip-40)
    else
      term2 = 0
    end

    z_bor_cb14=ztord_cb14 + w*Math.sin(dip*pi/180)

    zhypd_cb14 = magnitude + Math.exp([term1+term2,Math.log(0.9*(z_bor_cb14-magnitude))].min)
    if zhypd_cb14>20
      zhypd_cb14 = 0
    else
      zhypd_cb14 = [zhypd_cb14,5].max
    end

    wd_cb14 = 15; # these definitions are so roundabout. for the web site, the user must define w


    deltadpp = 0

    keyRegionLabel = {}
    keyRegionLabel[0] = "Global/California"
    keyRegionLabel[3] = "China"
    keyRegionLabel[4] = "Italy"
    keyRegionLabel[1] = "Japan"
    keyRegionLabel[2] = "Taiwan"
    keyRegionLabel[5] = "Turkey"

    regionLabel = keyRegionLabel[region]

    if regionLabel == "China"
      sj,sjp,sch,swn,sjpit,sji =   0,0,1,1,0,0
    elsif regionLabel == "Italy"
      sj,sjp,sch,swn,sjpit,sji =   1,1,0,0,1,1
    elsif regionLabel == "Japan"
      sj,sjp,sch,swn,sjpit,sji =   1,1,0,0,1,1
    elsif regionLabel == "Taiwan"
      sj,sjp,sch,swn,sjpit,sji =   0,0,0,0,0,0
    elsif regionLabel == "Turkey"
      sj,sjp,sch,swn,sjpit,sji =   0,0,1,1,0,0
    else        ## "California, Global"
      sj,sjp,sch,swn,sjpit,sji =   0,0,0,0,0,0
    end


    nGAModelSelection = params["NGAInputData_NGAModelSelection"]
    if (nGAModelSelection != "1")
      puts 'return'
      return
    end

    allGMPEstartCounter = 0

# -------------------------------------------------------------------------------------

    if(models[:AS] == true)
      thisGMPE = modelGMPEmapArray["AS"]
      # thisGMPE = "ask14"


      z1_to_ask14 = z10
      if z10 == 999
        if(region==1)
          set1 = (vs30*vs30+412*412)/(1360*1360+412*412)
          z1_to_ask14=Math.exp( -5.23/2*Math.log(set1))/1000
        else
          set2 = (vs30*vs30*vs30*vs30+610*610*610*610)/(1360*1360*1360*1360+610*610*610*610)
          z1_to_ask14=Math.exp( -7.67/4*Math.log(set2))/1000
        end
      end


      targetSpectrumGMPEArray["#{thisGMPE}median"]={}
      targetSpectrumGMPEArray["#{thisGMPE}stdev"]={}
      eval "thisWeight = w_#{thisGMPE}"
      targetSpectrumGMPEArray["#{thisGMPE}weight"]=thisWeight
      gmpeWeightSum += thisWeight
      gmpeList << thisGMPE


      thisGMPElabel = keyGMPELabelArray["#{thisGMPE}"]
      thisGMPEcolor = keyGMPEColorArray["#{thisGMPE}"]
      thisGMPEplotDataList = "{ lineWidth: 1, color: '#{thisGMPEcolor}', label: '#{thisGMPElabel}'}"
      if allGMPEstartCounter == 0
        allGMPEplotDataList = "#{thisGMPEplotDataList}"
      else
        allGMPEplotDataList = "#{allGMPEplotDataList},#{thisGMPEplotDataList}"
      end
      thisGMPEplotList = ""
      thisGMPEstartCounter = 0
      thisPeriodList = []
      thisMedianList = []
      thisSigmaList = []
      prepared_stmt = "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='rdgml_development' AND `TABLE_NAME`='z_gmpengaw2_#{thisGMPE}_coeffs' ORDER BY ORDINAL_POSITION; "
      result_setColumnName=(ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!
      counter = 0
      columnLabelList = []
      result_setColumnName.each{|thisColumnLabel|
        columnLabelList[counter] = thisColumnLabel.first.to_s
        counter+=1
      }
      prepared_stmt = "SELECT z_gmpengaw2_#{thisGMPE}_coeffs.* FROM z_gmpengaw2_#{thisGMPE}_coeffs;"
      result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!
      result_set.each do |coeffsRow|
        thisCounter = 0
        coeff = {}
        coeffsRow.each do |colValue|
          labelle = columnLabelList[thisCounter].to_s
          coeff[labelle] = colValue
          thisCounter +=1
        end
        unless coeff["period"]>0 then next end
        # ---------- the input in the following lines depends on the GMPE
        median = gmpe_ask_14(magnitude,rrup,rjb,rx,frv,fnm,fhw,fas,ztor,w,dip,vs30,vs30class,z1_to_ask14,ry0,region,coeff["period"],coeff["vlin"],coeff["b"],coeff["n"],coeff["m1"],coeff["c"],coeff["c4"],coeff["a1"],coeff["a2"],coeff["a3"],coeff["a4"],coeff["a5"],coeff["a6"],coeff["a7"],coeff["a8"],coeff["a10"],coeff["a11"],coeff["a12"],coeff["a13"],coeff["a14"],coeff["a15"],coeff["a16"],coeff["a17"],coeff["a43"],coeff["a44"],coeff["a45"],coeff["a46"],coeff["a25"],coeff["a28"],coeff["a29"],coeff["a31"],coeff["a36"],coeff["a37"],coeff["a38"],coeff["a39"],coeff["a40"],coeff["a41"],coeff["a42"],coeff["s1"],coeff["s2"],coeff["s3"],coeff["s4"],coeff["s1v"],coeff["s2v"],coeff["s5"],coeff["s6"])
        stdev = gmpe_ask14_stdev(magnitude,rrup,rjb,rx,frv,fnm,fhw,fas,ztor,w,dip,vs30,vs30class,z1_to_ask14,ry0,region,coeff["period"],coeff["vlin"],coeff["b"],coeff["n"],coeff["m1"],coeff["c"],coeff["c4"],coeff["a1"],coeff["a2"],coeff["a3"],coeff["a4"],coeff["a5"],coeff["a6"],coeff["a7"],coeff["a8"],coeff["a10"],coeff["a11"],coeff["a12"],coeff["a13"],coeff["a14"],coeff["a15"],coeff["a16"],coeff["a17"],coeff["a43"],coeff["a44"],coeff["a45"],coeff["a46"],coeff["a25"],coeff["a28"],coeff["a29"],coeff["a31"],coeff["a36"],coeff["a37"],coeff["a38"],coeff["a39"],coeff["a40"],coeff["a41"],coeff["a42"],coeff["s1"],coeff["s2"],coeff["s3"],coeff["s4"],coeff["s1v"],coeff["s2v"],coeff["s5"],coeff["s6"])

        # ---------- the input in the above lines depends on the GMPE
        targetSpectrumGMPEArray["#{thisGMPE}median"][coeff["period"]] = median
        targetSpectrumGMPEArray["#{thisGMPE}stdev"][coeff["period"]] = stdev
        thisPeriodList << coeff["period"]
        thisMedianList << median
        thisSigmaList << stdev
        if thisGMPEstartCounter == 0
          thisGMPEplotList = "[[#{coeff["period"]},#{median}]"
          thisGMPEstartCounter = 1
        else
          thisGMPEplotList = "#{thisGMPEplotList},[#{coeff["period"]},#{median}]"
        end
      end

      thisGMPEplotList = "#{thisGMPEplotList}]"
      if allGMPEstartCounter == 0
        allGMPEplotList = "#{thisGMPEplotList}"
        allGMPEstartCounter = 1
      else
        allGMPEplotList = "#{allGMPEplotList},#{thisGMPEplotList}"
      end
    end


# -------------------------------------------------------------------------------------
    if(models[:CB] == true)
      thisGMPE = modelGMPEmapArray["CB"]
      # thisGMPE = "cb14"



      targetSpectrumGMPEArray["#{thisGMPE}median"]={}
      targetSpectrumGMPEArray["#{thisGMPE}stdev"]={}
      eval "thisWeight = w_#{thisGMPE}"
      targetSpectrumGMPEArray["#{thisGMPE}weight"]=thisWeight
      gmpeWeightSum += thisWeight
      gmpeList << thisGMPE

      thisGMPElabel = keyGMPELabelArray["#{thisGMPE}"]
      thisGMPEcolor = keyGMPEColorArray["#{thisGMPE}"]
      thisGMPEplotDataList = "{ lineWidth: 1, color: '#{thisGMPEcolor}', label: '#{thisGMPElabel}'}"
      if allGMPEstartCounter == 0
        allGMPEplotDataList = "#{thisGMPEplotDataList}"
      else
        allGMPEplotDataList = "#{allGMPEplotDataList},#{thisGMPEplotDataList}"
      end
      thisGMPEplotList = ""
      thisGMPEstartCounter = 0
      thisPeriodList = []
      thisMedianList = []
      thisSigmaList = []
      prepared_stmt = "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='rdgml_development' AND `TABLE_NAME`='z_gmpengaw2_#{thisGMPE}_coeffs' ORDER BY ORDINAL_POSITION;"
      result_setColumnName=(ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!
      counter = 0
      columnLabelList = []
      result_setColumnName.each{|thisColumnLabel|
        columnLabelList[counter] = thisColumnLabel.first.to_s
        counter+=1
      }
      prepared_stmt = "SELECT z_gmpengaw2_#{thisGMPE}_coeffs.* FROM z_gmpengaw2_#{thisGMPE}_coeffs where period=0.001;"
      result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!
      a1100_cb = 0
      median0_cb = 0

      result_set.each do |coeffsRow|
        thisCounter = 0
        coeff = {}
        coeffsRow.each do |colValue|
          labelle = columnLabelList[thisCounter].to_s
          coeff[labelle] = colValue
          thisCounter +=1
        end
        # ---------- the input in the following lines depends on the GMPE
        a1100_cb = gmpe_a1100_cb(magnitude,rrup,rjb,rx,frv,fnm,fhw,w,dip,ztor,z25,zhyp,ztord_cb14,wd_cb14,zhypd_cb14,region,coeff["period"],coeff["c0"],coeff["c1"],coeff["c2"],coeff["c3"],coeff["c4"],coeff["c5"],coeff["c6"],coeff["c7"],coeff["c8"],coeff["c9"],coeff["c10"],coeff["c11"],coeff["c12"],coeff["c13"],coeff["c14"],coeff["c15"],coeff["c16"],coeff["c17"],coeff["c18"],coeff["c19"],coeff["c20"],coeff["dc20ca"],coeff["dc20jp"],coeff["dc20ch"],coeff["a2"],coeff["h1"],coeff["h2"],coeff["h3"],coeff["h4"],coeff["h5"],coeff["h6"],coeff["k1"],coeff["k2"],coeff["k3"],coeff["c"],coeff["n"],coeff["phi1"],coeff["phi2"],coeff["tau1"],coeff["tau2"],coeff["philnaf"],coeff["phic"],coeff["ro"],0)
        median0_cb = gmpe_cb_14(magnitude,rrup,rjb,rx,frv,fnm,fhw,ztor,w,dip,vs30,z25,zhyp,ztord_cb14,wd_cb14,zhypd_cb14,region,a1100_cb,coeff["period"],coeff["c0"],coeff["c1"],coeff["c2"],coeff["c3"],coeff["c4"],coeff["c5"],coeff["c6"],coeff["c7"],coeff["c8"],coeff["c9"],coeff["c10"],coeff["c11"],coeff["c12"],coeff["c13"],coeff["c14"],coeff["c15"],coeff["c16"],coeff["c17"],coeff["c18"],coeff["c19"],coeff["c20"],coeff["dc20ca"],coeff["dc20jp"],coeff["dc20ch"],coeff["a2"],coeff["h1"],coeff["h2"],coeff["h3"],coeff["h4"],coeff["h5"],coeff["h6"],coeff["k1"],coeff["k2"],coeff["k3"],coeff["c"],coeff["n"],coeff["phi1"],coeff["phi2"],coeff["tau1"],coeff["tau2"],coeff["philnaf"],coeff["phic"],coeff["ro"],0)
      end

      prepared_stmt = "SELECT z_gmpengaw2_#{thisGMPE}_coeffs.* FROM z_gmpengaw2_#{thisGMPE}_coeffs;"
      result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!
      result_set.each do |coeffsRow|
        thisCounter = 0
        coeff = {}
        coeffsRow.each do |colValue|
          labelle = columnLabelList[thisCounter].to_s
          coeff[labelle] = colValue
          thisCounter +=1
        end
        unless coeff["period"]>0 then next end
        # ---------- the input in the following lines depends on the GMPE
        median = gmpe_cb_14(magnitude,rrup,rjb,rx,frv,fnm,fhw,ztor,w,dip,vs30,z25,zhyp,ztord_cb14,wd_cb14,zhypd_cb14,region,a1100_cb,coeff["period"],coeff["c0"],coeff["c1"],coeff["c2"],coeff["c3"],coeff["c4"],coeff["c5"],coeff["c6"],coeff["c7"],coeff["c8"],coeff["c9"],coeff["c10"],coeff["c11"],coeff["c12"],coeff["c13"],coeff["c14"],coeff["c15"],coeff["c16"],coeff["c17"],coeff["c18"],coeff["c19"],coeff["c20"],coeff["dc20ca"],coeff["dc20jp"],coeff["dc20ch"],coeff["a2"],coeff["h1"],coeff["h2"],coeff["h3"],coeff["h4"],coeff["h5"],coeff["h6"],coeff["k1"],coeff["k2"],coeff["k3"],coeff["c"],coeff["n"],coeff["phi1"],coeff["phi2"],coeff["tau1"],coeff["tau2"],coeff["philnaf"],coeff["phic"],coeff["ro"],median0_cb)
        stdev = gmpe_cb_14stdev(magnitude,rrup,rjb,rx,frv,fnm,fhw,ztor,w,dip,vs30,z25,zhyp,ztord_cb14,wd_cb14,zhypd_cb14,region,a1100_cb,coeff["period"],coeff["c0"],coeff["c1"],coeff["c2"],coeff["c3"],coeff["c4"],coeff["c5"],coeff["c6"],coeff["c7"],coeff["c8"],coeff["c9"],coeff["c10"],coeff["c11"],coeff["c12"],coeff["c13"],coeff["c14"],coeff["c15"],coeff["c16"],coeff["c17"],coeff["c18"],coeff["c19"],coeff["c20"],coeff["dc20ca"],coeff["dc20jp"],coeff["dc20ch"],coeff["a2"],coeff["h1"],coeff["h2"],coeff["h3"],coeff["h4"],coeff["h5"],coeff["h6"],coeff["k1"],coeff["k2"],coeff["k3"],coeff["c"],coeff["n"],coeff["phi1"],coeff["phi2"],coeff["tau1"],coeff["tau2"],coeff["philnaf"],coeff["phic"],coeff["ro"])
        # ---------- the input in the above lines depends on the GMPE
        targetSpectrumGMPEArray["#{thisGMPE}median"][coeff["period"]] = median
        targetSpectrumGMPEArray["#{thisGMPE}stdev"][coeff["period"]] = stdev
        thisPeriodList << coeff["period"]
        thisMedianList << median
        thisSigmaList << stdev
        if thisGMPEstartCounter == 0
          thisGMPEplotList = "[[#{coeff["period"]},#{median}]"
          thisGMPEstartCounter = 1
        else
          thisGMPEplotList = "#{thisGMPEplotList},[#{coeff["period"]},#{median}]"
        end
      end
      thisGMPEplotList = "#{thisGMPEplotList}]"
      if allGMPEstartCounter == 0
        allGMPEplotList = "#{thisGMPEplotList}"
        allGMPEstartCounter = 1
      else
        allGMPEplotList = "#{allGMPEplotList},#{thisGMPEplotList}"
      end
    end

# -------------------------------------------------------------------------------------
    if(models[:BA] == true)
      # thisGMPE = "bssa14"
      thisGMPE = modelGMPEmapArray["BA"]


      targetSpectrumGMPEArray["#{thisGMPE}median"]={}
      targetSpectrumGMPEArray["#{thisGMPE}stdev"]={}
      eval "thisWeight = w_#{thisGMPE}"
      targetSpectrumGMPEArray["#{thisGMPE}weight"]=thisWeight
      gmpeWeightSum += thisWeight
      gmpeList << thisGMPE



      thisGMPElabel = keyGMPELabelArray["#{thisGMPE}"]
      thisGMPEcolor = keyGMPEColorArray["#{thisGMPE}"]
      thisGMPEplotDataList = "{ lineWidth: 1, color: '#{thisGMPEcolor}', label: '#{thisGMPElabel}'}"
      if allGMPEstartCounter == 0
        allGMPEplotDataList = "#{thisGMPEplotDataList}"
      else
        allGMPEplotDataList = "#{allGMPEplotDataList},#{thisGMPEplotDataList}"
      end
      thisGMPEplotList = ""
      thisGMPEstartCounter = 0
      thisPeriodList = []
      thisMedianList = []
      thisSigmaList = []

      #prepared_stmt = "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='rdgml_development' AND `TABLE_NAME`='z_gmpengaw2_#{thisGMPE}_coeffs';"

      prepared_stmt = "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='rdgml_development' AND `TABLE_NAME`='z_gmpengaw2_#{thisGMPE}_coeffs' ORDER BY ORDINAL_POSITION;"



      result_setColumnName=(ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!
      counter = 0
      columnLabelList = []
      result_setColumnName.each{|thisColumnLabel|
        columnLabelList[counter] = thisColumnLabel.first.to_s
        counter+=1
      }
      pgar_bssa = 0
      prepared_stmt = "SELECT z_gmpengaw2_#{thisGMPE}_coeffs.* FROM z_gmpengaw2_#{thisGMPE}_coeffs where period=0.0;"
      result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!
      result_set.each do |coeffsRow|
        thisCounter = 0
        coeff = {}
        coeffsRow.each do |colValue|
          labelle = columnLabelList[thisCounter].to_s
          coeff[labelle] = colValue
          thisCounter +=1
        end
        # ---------- the input in the following lines depends on the GMPE
        pgar_bssa = gmpe_pgar_calc(magnitude,rjb,u,frv,fnm,region,coeff["period"],coeff["e0"],coeff["e1"],coeff["e2"],coeff["e3"],coeff["e4"],coeff["e5"],coeff["e6"],coeff["mh"],coeff["c1"],coeff["c2"],coeff["c3"],coeff["mref"],coeff["rref"],coeff["h"],coeff["dc3"],coeff["dc3chtur"],coeff["dc3jpit"],coeff["c"],coeff["vc"],coeff["vref"],coeff["f1"],coeff["f3"],coeff["f4"],coeff["f5"],coeff["f6"],coeff["f7"],coeff["r1"],coeff["r2"],coeff["dfr"],coeff["dfv"],coeff["v1"],coeff["v2"],coeff["phi1"],coeff["phi2"],coeff["tau1"],coeff["tau2"])

      end
      prepared_stmt = "SELECT z_gmpengaw2_#{thisGMPE}_coeffs.* FROM z_gmpengaw2_#{thisGMPE}_coeffs;"
      result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!

      result_set.each do |coeffsRow|
        thisCounter = 0
        coeff = {}
        coeffsRow.each do |colValue|
          labelle = columnLabelList[thisCounter].to_s
          coeff[labelle] = colValue
          thisCounter +=1
        end
        unless coeff["period"]>0 then next end
        # ---------- the input in the following lines depends on the GMPE
        puts "pgar_bssa "+pgar_bssa.to_s
        puts "coeff[f3] "+coeff["f3"].to_s
        median = gmpe_bssa_14(magnitude,rjb,u,frv,fnm,vs30,region,z10,pgar_bssa,coeff["period"],coeff["e0"],coeff["e1"],coeff["e2"],coeff["e3"],coeff["e4"],coeff["e5"],coeff["e6"],coeff["mh"],coeff["c1"],coeff["c2"],coeff["c3"],coeff["mref"],coeff["rref"],coeff["h"],coeff["dc3"],coeff["dc3chtur"],coeff["dc3jpit"],coeff["c"],coeff["vc"],coeff["vref"],coeff["f1"],coeff["f3"],coeff["f4"],coeff["f5"],coeff["f6"],coeff["f7"],coeff["r1"],coeff["r2"],coeff["dfr"],coeff["dfv"],coeff["v1"],coeff["v2"],coeff["phi1"],coeff["phi2"],coeff["tau1"],coeff["tau2"])
        stdev = gmpe_bssa14_stdev(magnitude,rjb,vs30,coeff["period"],coeff["e0"],coeff["e1"],coeff["e2"],coeff["e3"],coeff["e4"],coeff["e5"],coeff["e6"],coeff["mh"],coeff["c1"],coeff["c2"],coeff["c3"],coeff["mref"],coeff["rref"],coeff["h"],coeff["dc3"],coeff["dc3chtur"],coeff["dc3jpit"],coeff["c"],coeff["vc"],coeff["vref"],coeff["f1"],coeff["f3"],coeff["f4"],coeff["f5"],coeff["f6"],coeff["f7"],coeff["r1"],coeff["r2"],coeff["dfr"],coeff["dfv"],coeff["v1"],coeff["v2"],coeff["phi1"],coeff["phi2"],coeff["tau1"],coeff["tau2"])
        # ---------- the input in the above lines depends on the GMPE
        targetSpectrumGMPEArray["#{thisGMPE}median"][coeff["period"]] = median
        targetSpectrumGMPEArray["#{thisGMPE}stdev"][coeff["period"]] = stdev
        thisPeriodList << coeff["period"]
        thisMedianList << median
        thisSigmaList << stdev

        if thisGMPEstartCounter == 0
          thisGMPEplotList = "[[#{coeff["period"]},#{median}]"
          thisGMPEstartCounter = 1
        else
          thisGMPEplotList = "#{thisGMPEplotList},[#{coeff["period"]},#{median}]"
        end
      end
      thisGMPEplotList = "#{thisGMPEplotList}]"
      if allGMPEstartCounter == 0
        allGMPEplotList = "#{thisGMPEplotList}"
        allGMPEstartCounter = 1
      else
        allGMPEplotList = "#{allGMPEplotList},#{thisGMPEplotList}"
      end
    end





# -------------------------------------------------------------------------------------
    if(models[:CY] == true)
      # thisGMPE = "cy14"
      thisGMPE = modelGMPEmapArray["CY"]



      z1r = 0
      if (region==1)
        z1r = Math.exp((-5.23/2)*Math.log((vs30**2+412**2)/(1360**2+412**2)))/1000
      else
        z1r = Math.exp((-7.15/4)*Math.log((vs30**4+571**4)/(1360**4+571**4)))/1000
      end
      ddpp = 0


      targetSpectrumGMPEArray["#{thisGMPE}median"]={}
      targetSpectrumGMPEArray["#{thisGMPE}stdev"]={}
      eval "thisWeight = w_#{thisGMPE}"
      targetSpectrumGMPEArray["#{thisGMPE}weight"]=thisWeight
      gmpeWeightSum += thisWeight
      gmpeList << thisGMPE



      thisGMPElabel = keyGMPELabelArray["#{thisGMPE}"]
      thisGMPEcolor = keyGMPEColorArray["#{thisGMPE}"]
      thisGMPEplotDataList = "{ lineWidth: 1, color: '#{thisGMPEcolor}', label: '#{thisGMPElabel}'}"
      if allGMPEstartCounter == 0
        allGMPEplotDataList = "#{thisGMPEplotDataList}"
      else
        allGMPEplotDataList = "#{allGMPEplotDataList},#{thisGMPEplotDataList}"
      end
      thisGMPEplotList = ""
      thisGMPEstartCounter = 0
      thisPeriodList = []
      thisMedianList = []
      thisSigmaList = []
      prepared_stmt = "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='rdgml_development' AND `TABLE_NAME`='z_gmpengaw2_#{thisGMPE}_coeffs' ORDER BY ORDINAL_POSITION;"
      result_setColumnName=(ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!
      counter = 0
      columnLabelList = []
      result_setColumnName.each{|thisColumnLabel|
        columnLabelList[counter] = thisColumnLabel.first.to_s
        counter+=1
      }
      ypga_cy = 0
      prepared_stmt = "SELECT z_gmpengaw2_#{thisGMPE}_coeffs.* FROM z_gmpengaw2_#{thisGMPE}_coeffs where period=0;"
      result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!
      result_set.each do |coeffsRow|
        thisCounter = 0
        coeff = {}
        coeffsRow.each do |colValue|
          labelle = columnLabelList[thisCounter].to_s
          coeff[labelle] = colValue
          # statm = "#{labelle} = #{colValue}"
          # eval(statm)
          thisCounter +=1
        end
        # ---------- the input in the following lines depends on the GMPE
        # v["ypga_cy"] = gmpe_cy_14(v["m"],v["rrup"],v["rjb"],v["rx"],v["vs30"],v["frv"],v["fnm"],v["fhw"],v["dip"],v["ztor_user"],v["region"],v["z10user"],v["z1r"],v["ddpp"],coeff["period"],coeff["c2"],coeff["c4"],coeff["c4a"],coeff["crb"],coeff["c8"],coeff["c8a"],coeff["c1"],coeff["c1a"],coeff["c1b"],coeff["c1c"],coeff["c1d"],coeff["cn"],coeff["cm"],coeff["c3"],coeff["c5"],coeff["chm"],coeff["c6"],coeff["c7"],coeff["c7b"],coeff["c8b"],coeff["c9"],coeff["c9a"],coeff["c9b"],coeff["c11"],coeff["c11b"],coeff["cg1"],coeff["cg2"],coeff["cg3"],coeff["f1"],coeff["f2"],coeff["f3"],coeff["f4"],coeff["f5"],coeff["f6"],coeff["t1"],coeff["t2"],coeff["s1"],coeff["s2"],coeff["s3"],coeff["s2jp"],coeff["gjpit"],coeff["gwn"],coeff["f1jp"],coeff["f5jp"],coeff["f6jp"],0)
        ypga_cy = gmpe_cy_14(magnitude,rrup,rjb,rx,vs30,frv,fnm,fhw,dip,ztor,region,z10,z1r,ddpp,coeff["period"],coeff["c2"],coeff["c4"],coeff["c4a"],coeff["crb"],coeff["c8"],coeff["c8a"],coeff["c1"],coeff["c1a"],coeff["c1b"],coeff["c1c"],coeff["c1d"],coeff["cn"],coeff["cm"],coeff["c3"],coeff["c5"],coeff["chm"],coeff["c6"],coeff["c7"],coeff["c7b"],coeff["c8b"],coeff["c9"],coeff["c9a"],coeff["c9b"],coeff["c11"],coeff["c11b"],coeff["cg1"],coeff["cg2"],coeff["cg3"],coeff["f1"],coeff["f2"],coeff["f3"],coeff["f4"],coeff["f5"],coeff["f6"],coeff["t1"],coeff["t2"],coeff["s1"],coeff["s2"],coeff["s3"],coeff["s2jp"],coeff["gjpit"],coeff["gwn"],coeff["f1jp"],coeff["f5jp"],coeff["f6jp"],0)
      end
      prepared_stmt = "SELECT z_gmpengaw2_#{thisGMPE}_coeffs.* FROM z_gmpengaw2_#{thisGMPE}_coeffs;"
      result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!
      result_set.each do |coeffsRow|
        thisCounter = 0
        coeff = {}
        coeffsRow.each do |colValue|
          labelle = columnLabelList[thisCounter].to_s
          coeff[labelle] = colValue
          thisCounter +=1
        end
        unless coeff["period"]>0 then next end
        # ---------- the input in the following lines depends on the GMPE
        median = gmpe_cy_14(magnitude,rrup,rjb,rx,vs30,frv,fnm,fhw,dip,ztor,region,z10,z1r,ddpp,coeff["period"],coeff["c2"],coeff["c4"],coeff["c4a"],coeff["crb"],coeff["c8"],coeff["c8a"],coeff["c1"],coeff["c1a"],coeff["c1b"],coeff["c1c"],coeff["c1d"],coeff["cn"],coeff["cm"],coeff["c3"],coeff["c5"],coeff["chm"],coeff["c6"],coeff["c7"],coeff["c7b"],coeff["c8b"],coeff["c9"],coeff["c9a"],coeff["c9b"],coeff["c11"],coeff["c11b"],coeff["cg1"],coeff["cg2"],coeff["cg3"],coeff["f1"],coeff["f2"],coeff["f3"],coeff["f4"],coeff["f5"],coeff["f6"],coeff["t1"],coeff["t2"],coeff["s1"],coeff["s2"],coeff["s3"],coeff["s2jp"],coeff["gjpit"],coeff["gwn"],coeff["f1jp"],coeff["f5jp"],coeff["f6jp"],ypga_cy)
        stdev = gmpe_cy14_stdev(magnitude,rrup,rjb,rx,vs30,frv,fnm,fhw,dip,ztor,region,z10,ddpp,vs30class,coeff["period"],coeff["c2"],coeff["c4"],coeff["c4a"],coeff["crb"],coeff["c8"],coeff["c8a"],coeff["c1"],coeff["c1a"],coeff["c1b"],coeff["c1c"],coeff["c1d"],coeff["cn"],coeff["cm"],coeff["c3"],coeff["c5"],coeff["chm"],coeff["c6"],coeff["c7"],coeff["c7b"],coeff["c8b"],coeff["c9"],coeff["c9a"],coeff["c9b"],coeff["c11"],coeff["c11b"],coeff["cg1"],coeff["cg2"],coeff["cg3"],coeff["f1"],coeff["f2"],coeff["f3"],coeff["f4"],coeff["f5"],coeff["f6"],coeff["t1"],coeff["t2"],coeff["s1"],coeff["s2"],coeff["s3"],coeff["s2jp"],coeff["gjpit"],coeff["gwn"],coeff["f1jp"],coeff["f5jp"],coeff["f6jp"])
        # ---------- the input in the above lines depends on the GMPE
        targetSpectrumGMPEArray["#{thisGMPE}median"][coeff["period"]] = median
        targetSpectrumGMPEArray["#{thisGMPE}stdev"][coeff["period"]] = stdev
        thisPeriodList << coeff["period"]
        thisMedianList << median
        thisSigmaList << stdev
        if thisGMPEstartCounter == 0
          thisGMPEplotList = "[[#{coeff["period"]},#{median}]"
          thisGMPEstartCounter = 1
        else
          thisGMPEplotList = "#{thisGMPEplotList},[#{coeff["period"]},#{median}]"
        end
      end

      thisGMPEplotList = "#{thisGMPEplotList}]"
      if allGMPEstartCounter == 0
        allGMPEplotList = "#{thisGMPEplotList}"
        allGMPEstartCounter = 1
      else
        allGMPEplotList = "#{allGMPEplotList},#{thisGMPEplotList}"
      end
    end
# -------------------------------------------------------------------------------------
    if(models[:ID] == true)
      # thisGMPE = "i14"
      thisGMPE = modelGMPEmapArray["ID"]

      targetSpectrumGMPEArray["#{thisGMPE}median"]={}
      targetSpectrumGMPEArray["#{thisGMPE}stdev"]={}
      eval "thisWeight = w_#{thisGMPE}"
      targetSpectrumGMPEArray["#{thisGMPE}weight"]=thisWeight
      gmpeWeightSum += thisWeight
      gmpeList << thisGMPE



      thisGMPElabel = keyGMPELabelArray["#{thisGMPE}"]
      thisGMPEcolor = keyGMPEColorArray["#{thisGMPE}"]
      thisGMPEplotDataList = "{ lineWidth: 1, color: '#{thisGMPEcolor}', label: '#{thisGMPElabel}'}"
      if allGMPEstartCounter == 0
        allGMPEplotDataList = "#{thisGMPEplotDataList}"
      else
        allGMPEplotDataList = "#{allGMPEplotDataList},#{thisGMPEplotDataList}"
      end
      thisGMPEplotList = ""
      thisGMPEstartCounter = 0
      thisPeriodList = []
      thisMedianList = []
      thisSigmaList = []
      prepared_stmt = "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='rdgml_development' AND `TABLE_NAME`='z_gmpengaw2_#{thisGMPE}_coeffs' ORDER BY ORDINAL_POSITION;"
      result_setColumnName=(ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!
      counter = 0
      columnLabelList = []
      result_setColumnName.each{|thisColumnLabel|
        columnLabelList[counter] = thisColumnLabel.first.to_s
        counter+=1
      }
      prepared_stmt = "SELECT z_gmpengaw2_#{thisGMPE}_coeffs.* FROM z_gmpengaw2_#{thisGMPE}_coeffs;"
      result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!
      result_set.each do |coeffsRow|
        thisCounter = 0
        coeff = {}
        coeffsRow.each do |colValue|
          labelle = columnLabelList[thisCounter].to_s
          coeff[labelle] = colValue
          thisCounter +=1
        end
        unless coeff["period"]>0 then next end
        # ---------- the input in the following lines depends on the GMPE
        median = gmpe_i_14(magnitude,rrup,frv,vs30,coeff["period"],coeff["a11"],coeff["a21"],coeff["b11"],coeff["b21"],coeff["a12"],coeff["a22"],coeff["b12"],coeff["b22"],coeff["a3"],coeff["ksi"],coeff["g"],coeff["phi"])
        stdev = gmpe_i_14_stdev(magnitude,coeff["period"],coeff["a11"],coeff["a21"],coeff["b11"],coeff["b21"],coeff["a12"],coeff["a22"],coeff["b12"],coeff["b22"],coeff["a3"],coeff["ksi"],coeff["g"],coeff["phi"])
        # ---------- the input in the above lines depends on the GMPE
        targetSpectrumGMPEArray["#{thisGMPE}median"][coeff["period"]] = median
        targetSpectrumGMPEArray["#{thisGMPE}stdev"][coeff["period"]] = stdev
        thisPeriodList << coeff["period"]
        thisMedianList << median
        thisSigmaList << stdev
        if thisGMPEstartCounter == 0
          thisGMPEplotList = "[[#{coeff["period"]},#{median}]"
          thisGMPEstartCounter = 1
        else
          thisGMPEplotList = "#{thisGMPEplotList},[#{coeff["period"]},#{median}]"
        end
      end
      thisGMPEplotList = "#{thisGMPEplotList}]"
      if allGMPEstartCounter == 0
        allGMPEplotList = "#{thisGMPEplotList}"
        allGMPEstartCounter = 1
      else
        allGMPEplotList = "#{allGMPEplotList},#{thisGMPEplotList}"
      end
    end

# -----------------------------------------





    sql3 = "SELECT z_gmpe2_ask13coeffs.period
              FROM z_gmpe2_ask13coeffs INNER JOIN z_gmpe2_bssa13coeffs ON z_gmpe2_ask13coeffs.period = z_gmpe2_bssa13coeffs.period
               INNER JOIN z_gmpe2_i13coeffs ON z_gmpe2_bssa13coeffs.period = z_gmpe2_i13coeffs.period
               INNER JOIN z_gmpe2_cy13coeffs ON z_gmpe2_i13coeffs.period = z_gmpe2_cy13coeffs.period
               INNER JOIN z_gmpe2_cb13coeffs ON z_gmpe2_cy13coeffs.period = z_gmpe2_cb13coeffs.period
               AND z_gmpe2_cb13coeffs.period = z_gmpe2_ask13coeffs.period AND z_gmpe2_ask13coeffs.period = z_gmpe2_cy13coeffs.period
               AND z_gmpe2_ask13coeffs.period = z_gmpe2_i13coeffs.period AND z_gmpe2_bssa13coeffs.period = z_gmpe2_cy13coeffs.period
               AND z_gmpe2_bssa13coeffs.period = z_gmpe2_cb13coeffs.period AND z_gmpe2_i13coeffs.period = z_gmpe2_cb13coeffs.period"

    result_set=(ActiveRecord::Base.connection.execute(sql3))
    ActiveRecord::Base.connection_handler.clear_active_connections!

    gmpePeriodList = []
    epsilonArray = {}
    result_set.each{|thisRow|
      thisT = thisRow[0]
      next unless thisT>=0.01
      gmpePeriodList << thisT
      epsilonArray[thisT] = 1
    }

    if cmflag == 1
      epsilonArray = compute_epsilonArray(tdesign,gmpePeriodList)
    end

    targetpSaList = []
    targetPeriodList = []
    targetpSaArray = {}
    gmpePeriodList.each{|thisT|
      unless thisT>0 then next end

      thispSa = 0
      thisLogpSa = 0
      gmpeList.each{|thisGMPE|
        thisMedian = targetSpectrumGMPEArray["#{thisGMPE}median"][thisT]
        thisSigma = targetSpectrumGMPEArray["#{thisGMPE}stdev"][thisT]
        thisWeight = targetSpectrumGMPEArray["#{thisGMPE}weight"]/gmpeWeightSum
        thisEpsilon = epsilonArray[thisT]
        thisMedianPlusEpsilonSigma = Math.exp(Math.log(thisMedian) + epstdesign*thisEpsilon*thisSigma)
        thisWeightedMedianPlusEpsilonSigma = thisWeight*(thisMedianPlusEpsilonSigma)
        thispSa += thisWeightedMedianPlusEpsilonSigma
        thisLogpSa += thisWeight*(Math.log(thisMedianPlusEpsilonSigma))
      };
      thisLogpSa = Math.exp(thisLogpSa)


      if ngameantype == "Geometric"
        targetpSa = thisLogpSa
      else
        targetpSa = thispSa
      end
      targetPeriodList << thisT
      targetpSaList << targetpSa
      targetpSaArray[thisT] = targetpSa

    }

    targetSpectrumGMPEArray["targetPeriodList"] = targetPeriodList
    targetSpectrumGMPEArray["targetpSaList"] = targetpSaList
    targetSpectrumGMPEArray["targetpSaArray"] = targetpSaArray
    targetSpectrumGMPEArray["allGMPEplotList"] = allGMPEplotList
    targetSpectrumGMPEArray["allGMPEplotDataList"] = allGMPEplotDataList

    return targetSpectrumGMPEArray

  end


  def spectraImage
    SpectraImage.new(self)
  end


  def getRandomString
    #Time.now.to_f.to_s
    self.id.to_s
  end

  def convertStringToArray(str)
    brackets_removed_str = str[1,str.length-2]
    splitToArray = brackets_removed_str.split(',')

    array = []

    #convert the strings ["1", "1"...] to floats, caused error without this
    splitToArray.each do | a |
      array << a.to_f
    end

    return array
  end

  def validate_fields_user_defined_spectrum?()
    if self.NGAInputData_NGAModelSelection != 0 && self.NGAInputData_NGAModelSelection != 99 && self.NGAInputData_NGAModelSelection != 88
      return false;
    end
  end

  def validate_fields_asce_sds?(ngaModelSelection, s_sds)
    if ngaModelSelection == 99
      if s_sds.to_s.empty?
        return false
      end
    end
    return true
  end
  def validate_fields_asce_sd1?(ngaModelSelection,s_sd1)
    if ngaModelSelection == 99
      if s_sd1.to_s.empty?
        return false
      end
    end
    return true
  end
  def validate_fields_asce_TL?(ngaModelSelection, s_TL)
    if ngaModelSelection == 99
      if s_TL.to_s.empty?
        return false
      end
    end
    return true
  end

  def validate_fields_vs30?(as, ba, cb, cy, id, input_str2, vs30)
    #1. get checkbox String [1,1,1,1,1]
    #2. get list of valid models, make sure one model is checked
    #3. compare list against checkbox string for each valid model
    #4. if any one of the models are checked, then we must validate
    #   and return true from here

    if input_str2 == nil || vs30 == nil
      return false;
    end

    if self.NGAInputData_NGAModelSelection == 0
      return false;
    end

    brackets_removed_str = input_str2[1,input_str2.length-2]
    array = brackets_removed_str.split(',')


    if(array[4] == id.to_s && id == 1)
      if(vs30 >= 450)
        return true;
      else
        errors.add("Error", "Idriss is NOT applicable to cases where Vs30 < 450 m/s or Normal Faulting");
      end

      return false;
    end

    return false;
  end



  def validate_field_gmpe?(input_str2)
    if (self.NGAInputData_NGAModelSelection == 1 )
      if input_str2 == nil
        return false;
      end
    end

  end
  def validate_fields?(as, ba, cb, cy, id, input_str2)
    #1. get checkbox String [1,1,1,1,1]
    #2. get list of valid models, make sure one model is checked
    #3. compare list against checkbox string for each valid model
    #4. if any one of the models are checked, then we must validate
    #   and return true from here

    if input_str2 == nil
      return false;
    end

    if self.NGAInputData_NGAModelSelection == 0
      return false;
    end



    brackets_removed_str = input_str2[1,input_str2.length-2]
    array = brackets_removed_str.split(',')

    #check AS
    if(array[0] == as.to_s && as == 1)
      return true;
    end

    if(array[1] == ba.to_s && ba == 1)
      return true;
    end

    if(array[2] == cb.to_s && cb == 1)
      return true;
    end

    if(array[3] == cy.to_s && cy == 1)
      return true;
    end

    if(array[4] == id.to_s && id == 1)
      return true;
    end

    return false;
  end

  def validate_fields_Tdesign?(as, ba, cb, cy, id, cmflag, designT)

    if self.NGAInputData_NGAModelSelection == 0
      return false;
    end

    if cmflag == 1
      if designT == nil
        errors.add("Error", "The Conditioning Period for the CMS must be defined");
        return false;
      end
      if designT >=0.01 && designT <= 10
        return true;
      else
        errors.add("Error", "The Conditioning Period for the CMS must be between 0.01 and 10");
        return false
      end
    end

    return false;
  end

  def validate_fields_EpsTdesign?(as, ba, cb, cy, id, cmflag)

    if self.NGAInputData_NGAModelSelection == 0
      return false;
    end

    if cmflag == 1
      return true
    end

    return false;
  end
  def validate
    if (self.NGAInputData_checkboxNGA.include?("1") == false &&
        self.NGAInputData_NGAModelSelection != 0 &&
        self.NGAInputData_NGAModelSelection != 99 &&
        self.NGAInputData_NGAModelSelection != 88 )
      errors.add("Please", " select at least one model")
    end

    if(self.NGAInputData_NGAModelSelection == 0 && self.filename == "")
      errors.add("Please", " enter User Defined Spectrum filename")
    end

    if(self.NGAInputData_NGAModelSelection == 99 &&
        (self.NGAInputData_Sds.to_s.empty? ||
            self.NGAInputData_Sd1.to_s.empty? ||
            self.NGAInputData_TL.to_s.empty?)
    )
      errors.add("Please", " input ASCE Code Specification values")
    end


  end


  def self.human_attribute_name(attribute)
    #HUMANIZED_COLLUMNS[attribute.to_sym] || super
  end




  def CreateSpectraImage_WriteSpectraResults(params,targetSpectrumArray)

    puts "***********************"
    puts "***********************"
    puts " "
    puts 'Start in spectras model ->  CreateSpectraImage_WriteSpectraResults(params,targetSpectrumArray)'
    puts " "
    puts "***********************"
    puts "***********************"


    # startTime=Time.now
    unique_token=self.id.to_s

    thisSpectraID = self.id


    outputLines = []



    inList = [
        "",
        "  -- PEER Ground Motion Database Target Spectrum Report -- NGA-West2 -- " + Date.today.to_s,
        "",
    ]
    inList.each do |prepared_stmt|
      ##x## fout.puts(prepared_stmt)
      outputLines << prepared_stmt
    end


    # "No Target Spectrum"-88,"PEER-NGA West2 Spectrum"-1, "User Defined Spectrum"-0, "ASCE Code Spectrum"-99

    thisNGAModelSelection = params[:NGAInputData_NGAModelSelection]



    if ( thisNGAModelSelection == 88 || thisNGAModelSelection == "88")


      outputLines << ' -- Target Spectrum Type:  No Target Spectrum --'
      calcTime = Time.now-startTime
      puts 'writeSpectraResults - processing time: ' + calcTime.to_s

      return 1
    elsif (thisNGAModelSelection == 99 || thisNGAModelSelection == "99")


      outputLines << ' -- Target Spectrum Type:  ASCE Code Spectrum --'



      outputLines << ' -- Spectrum Parameters: --'
      outputLines << 'Sds (g):,' + self.NGAInputData_Sds.to_s
      outputLines << 'Sd1 (g):,' + self.NGAInputData_Sd1.to_s
      outputLines << 'TL (sec):,' + self.NGAInputData_TL.to_s

    elsif (thisNGAModelSelection == 0 || thisNGAModelSelection == "0")
      outputLines << ' -- Target Spectrum Type:  User Defined Spectrum --'
      fullFilename = params[:filename]
      puts 'fullFilename'
      puts fullFilename
      endIndex = fullFilename.index('#x#')
      userFilename = fullFilename[0..endIndex-1]
      outputLines << ' Input Filename:  '+userFilename.to_s
    else

      outputLines << ' -- Target Spectrum Type:  PEER NGA-West2 Spectrum --'
      inList = [
          "",
          " -- Earthquake Source Parameters for GMPE: -- ",
      ]
      inList.each do |prepared_stmt|

        outputLines << prepared_stmt
      end




      prepared_stmt = "SELECT id,NGAInputData_NGAModelSelection,NGAInputData_checkboxNGA,NGAInputData_Nsigma,menu_Mechanism,NGAInputData_M,NGAInputData_RRUP,NGAInputData_RJB,NGAInputData_Rx,NGAInputData_HWflag,NGAInputData_ZTOR,NGAInputData_ZHYP,NGAInputData_DIP,NGAInputData_Vs30,NGAInputData_Vs30_class,NGAInputData_Z10,NGAInputData_Z25,NGAInputData_W,NGAInputData_RY0,NGAInputData_Region,NGAInputData_CMFlag,NGAInputData_EpsTdesign,NGAInputData_Tdesign,temp,NGAMeanFlag,varargin,unique_token,xyScale,filename,NGAInputData_DampingRatio,w_ask14,w_bssa14,w_cb14,w_cy14,w_i14
                        FROM spectras where (id = #{thisSpectraID});"
      result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!

      id,ngamodelselection,checkboxnga,nsigma,menu_mechanism,magnitude,rrup,rjb,rx,hwflag,ztor,zhyp,dip,vs30,vs30_class,z10,z25,w,ry0,region,cmflag,epstdesign,tdesign,temp,ngameanflag,varargin,unique_token,xyscale,filename,dampingratio,w_ask14,w_bssa14,w_cb14,w_cy14,w_i14 = result_set.first
      keyNGAMeanArray={}
      keyNGAMeanArray[0] = "Geometric"
      keyNGAMeanArray[1] = "Arithmetic"
      ngameantype = keyNGAMeanArray[ngameanflag]


      keyMenuMechanismArray={}
      keyMenuMechanismArray[1]="Strike Slip"
      keyMenuMechanismArray[2]="Normal/Oblique"
      keyMenuMechanismArray[3]="Reverse/Oblique"
      faulttype = keyMenuMechanismArray[menu_mechanism.to_i]


      checkboxngaList = checkboxnga[1..-2].split(",")  # the range removes the brackets on the outside

      vs30type = "Estimated"
      if vs30_class == 0
        vs30type = "Measured"
      end
      hangingwalltype = "No"
      if hwflag == 1
        hangingwalltype = "Yes"
      end
      cmstype = "No"
      if cmflag == 1
        cmstype = "Yes"
      end
      theseLines = ",Abrahamson-Silva-Kamai '14:,#{checkboxngaList[0]}
                    ,Boore-Stewart-Seyhan-Atkinson '14:,#{checkboxngaList[1]}
                    ,Campbell-Bozorgnia '14:,#{checkboxngaList[2]}
                    ,Chiou-Youngs '14:,#{checkboxngaList[3]}
                    ,Idriss '14:,#{checkboxngaList[4]}
                    ,Weight - Abrahamson-Silva-Kamai '14:,#{w_ask14}
                    ,Weight - Boore-Stewart-Seyhan-Atkinson '14:,#{w_bssa14}
                    ,Weight - Campbell-Bozorgnia '14:,#{w_cb14}
                    ,Weight - Chiou-Youngs '14:,#{w_cy14}
                    ,Weight - Idriss '14:,#{w_i14}
                    ,Damping Ratio:,#{dampingratio}
                    ,Region:,#{region}
                    ,Fault Type:,#{faulttype}
                    ,Magnitude:,#{magnitude}
                    ,Rrup (km):,#{rrup}
                    ,Rx (km):,#{rx}
                    ,Ry0 (km):,#{ry0}
                    ,Rjb (km):,#{rjb}
                    ,Ztor (km):,#{ztor}
                    ,Width (km):,#{w}
                    ,DIP (deg):,#{dip}
                    ,Vs30 (m/s):,#{vs30}
                    ,Vs30 type:,#{vs30type}
                    ,Z1.0 (km):,#{z10}
                    ,Z2.5 (km):,#{z25}
                    ,Zhyp (km):,#{zhyp}
                    ,Hanging Wall:,#{hangingwalltype}
                    ,GMPE Average:,#{ngameantype}
                    ,Epsilon:,#{epstdesign}
                    ,CMS:,#{cmstype}
                    ,Tcms(s):,#{tdesign}
      "


      outputLines << theseLines


      outputLines << ""
      outputLines << ' -- GMPE Median Spectra --'

      thisRowString = "Period (sec)"
      if defined? targetSpectrumArray["ask14median"][0.01]
        thisRowString = "#{thisRowString},Abrahamson-Silva-Kamai '14"
      end
      if defined? targetSpectrumArray["bssa14median"][0.01]
        thisRowString = "#{thisRowString},Boore-Stewart-Seyhan-Atkinson '14"
      end
      if defined? targetSpectrumArray["cb14median"][0.01]
        thisRowString = "#{thisRowString},Campbell-Bozorgnia '14"
      end
      if defined? targetSpectrumArray["cy14median"][0.01]
        thisRowString = "#{thisRowString},Chiou-Youngs '14 '14"
      end
      if defined? targetSpectrumArray["i14median"][0.01]
        thisRowString = "#{thisRowString},Idriss '14"
      end
      outputLines << thisRowString

      targetPeriodList = targetSpectrumArray["targetPeriodList"]
      targetpSaArray = targetSpectrumArray["targetpSaArray"]




      targetPeriodList.each {|thisT|
        thisRowString = "#{thisT}"
        if defined? targetSpectrumArray["ask14median"][thisT]
          thisRowString = "#{thisRowString},#{targetSpectrumArray["ask14median"][thisT].round(8)}"
        end
        if defined? targetSpectrumArray["bssa14median"][thisT]
          thisRowString = "#{thisRowString},#{targetSpectrumArray["bssa14median"][thisT].round(8)}"
        end
        if defined? targetSpectrumArray["cb14median"][thisT]
          thisRowString = "#{thisRowString},#{targetSpectrumArray["cb14median"][thisT].round(8)}"
        end
        if defined? targetSpectrumArray["cy14median"][thisT]
          thisRowString = "#{thisRowString},#{targetSpectrumArray["cy14median"][thisT].round(8)}"
        end
        if defined? targetSpectrumArray["i14median"][thisT]
          thisRowString = "#{thisRowString},#{targetSpectrumArray["i14median"][thisT].round(8)}"
        end
        outputLines << thisRowString

      }

    end


    outputLines << ""
    outputLines << ' -- Target Spectrum --'
    outputLines << 'T (sec),pSa (g)'



    targetPeriodList = targetSpectrumArray["targetPeriodList"]
    targetpSaList = targetSpectrumArray["targetpSaList"]
    targetpSaArray = targetSpectrumArray["targetpSaArray"]

    targetPeriodList.each{|thisT|
      thispSa = targetpSaArray[thisT]
      outputLines << thisT.to_s + ',' + thispSa.round(8).to_s
    }


    thisFilename = "TargetSpectraResults#{unique_token}.csv"


    # directory = "/archive/ngawest2_files/Utility"
    # directory = "#{ENV['ARCHIVE_FILES_PATH']}Utility"
    directory_download = "#{ENV['ARCHIVE_FILES_DOWNLOAD_PATH']}Utility"

    thisFullPath = File.join(directory_download, thisFilename)
    thisFullPath = File.expand_path(thisFullPath)



    fout = File.open(thisFullPath, 'w')

    fout.puts(outputLines)

    fout.close




    # calcTime = Time.now-startTime


    return 1





  end







  def compute_epsilonArray(tDesign, tList)
    # given array of tList, and a single design period tDesign, find the
    # correlation array corresponding to each t
    # user baker's equation
    pi = 4 * Math.atan(1)
    epsilonArray = {}
    tList.each{|thisT|
      if thisT> tDesign
        tmin=tDesign;
        tmax=thisT;
      else
        tmin=thisT;
        tmax=tDesign;
      end
      ##****************************************************************
      ## use baker & jayaram: correlation of spectral acceleration values
      ## from nga ground motion models, eq spectra, 2008
      ##****************************************************************
      c1=1-Math.cos(pi/2-0.366*Math.log(tmax/[tmin, 0.109].max));
      if tmax<0.2
        c2=1-0.105*(1.0-1.0/(1.0+Math.exp(100*tmax-5)))*(tmax-tmin)/(tmax-0.0099);
      else
        c2=0;
      end
      if tmax<0.109
        c3=c2;
      else
        c3=c1;
      end
      c4=c1+0.5*(Math.sqrt(c3)-c3)*(1+Math.cos(pi*tmin/0.109));
      if tmax<0.109            # for lower left corner
        rho_xx=c2;
      elsif tmin>0.109        # for upper right corner
        rho_xx=c1;
      elsif tmax<0.2
        rho_xx=[c2, c4].min;
      else
        rho_xx=c4;
      end
      epsilonArray[thisT] = rho_xx
    }
    return epsilonArray
  end

  ##


  #######################################################
  #######################################################
  ################ NGA-WEST2 GMPEs#######################
  #######################################################
  #######################################################

  ## public ypgacb as double

  def gmpe_a1100_cb(m,rrup,rjb,rx,frv,fnm,fhw,w,del,ztor,z25,zhyp,ztord,wd,zhypd,region,t,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,dc20ca,dc20jp,dc20ch,a2,h1,h2,h3,h4,h5,h6,k1,k2,k3,c,n,phi1,phi2,tau1,tau2,philnaf,phic,ro,ypgacb)
    # ***************************************************************
    # *                                                             *
    # *      campbell-bozorgnia nga model for pgar, eqs versio      *
    # *               written by emel seyhan, phd                   *
    # *                  last checked 03 23 14                      *
    # *                                                             *
    # ***************************************************************
    # use coeffs for pga
    ## dim counter as integer
    # counter = 0
    # ## dim cof(36) as double
    # case t
    # when = 0
    #     for each cell in range("cb14_coeffs!b26:cb14_coeffs!ak26")
    #         cof(counter) = cell.value
    #         counter = counter + 1
    #     next cell
    # end
    # c0 = cof(0)
    # c1 = cof(1)
    # c2 = cof(2)
    # c3 = cof(3)
    # c4 = cof(4)
    # c5 = cof(5)
    # c6 = cof(6)
    # c7 = cof(7)
    # n = cof(35)
    pi = 4 * Math.atan(1)
    if (region == 1)
      sj = 1
    else
      sj = 0
    end

    if (ztor == 999 || ztor == nil)
      ztor = ztord
    else
      ztor = ztor
    end
    # default hypocentral depth
    if (zhyp == 999)
      zhyp = zhypd
    else
      zhyp = zhyp
    end
    if (w == 999)
      zhyp = 9
    else
      zhyp = zhyp
    end
    # default w
    if (w == 999)
      w = wd
    else
      w = w
    end
    # reference z2.5 # ca and jp
    z25r = (1 - sj) * Math.exp(7.089 - 1.144 * Math.log(1100)) + sj * Math.exp(5.359 - 1.102 * Math.log(1100))
    if (z25 == 999)
      z25 = z25r
    else
      z25 = z25
    end
    # magnitude scaling term
    case
      when m <= 4.5
        fmag = (c1 * m)
      when m <= 5.5
        fmag = (c1 * m) + (c2 * (m - 4.5))
      when m <= 6.5
        fmag = (c1 * m) + (c2 * (m - 4.5)) + (c3 * (m - 5.5))
      else
        fmag = (c1 * m) + (c2 * (m - 4.5)) + (c3 * (m - 5.5)) + (c4 * (m - 6.5))
    end
    # distance scaling term
    fdis = ((c5 + (c6 * m)) * (Math.log(Math.sqrt((rrup * rrup) + (c7 * c7)))))
    # style of faulting term
    ffltf = (c8 * frv) + (c9 * fnm)
    case
      when m <= 4.5
        ffltm = 0
      when m <= 5.5
        ffltm = m - 4.5
      else
        ffltm = 1
    end
    fflt = ffltf * ffltm
    # hanging wall term

    r1 = w * Math.cos(del * pi / 180)
    r2 = 62 * m - 350
    f1_rx = h1 + h2 * (rx / r1) + h3 * (rx / r1)**2
    f2_rx = h4 + h5 * ((rx - r1) / (r2 - r1)) + h6 * ((rx - r1) / (r2 - r1))**2
    if (f2_rx > 0)
      maxr = f2_rx
    else
      maxr = 0
    end
    if (fhw == 0)
      fhngrx = 0
    elsif (rx < r1 and fhw == 1)
      fhngrx = f1_rx
    elsif (rx >= r1 and fhw == 1)
      fhngrx = maxr
    end
    case
      when rrup <= 0
        fhngr = 1
      else
        fhngr = (rrup - rjb) / rrup
    end
    if (f2_rx > 0)
      f2x = f2_rx
    else
      f2x = 0
    end
    case
      when m <= 5.5
        fhngm = 0
      when m <= 6.5
        fhngm = (m - 5.5) * (1 + a2 * (m - 6.5))
      else
        fhngm = 1 + a2 * (m - 6.5)
    end
    case
      when ztor < 16.66
        fhngz = 1 - 0.06 * ztor
      else
        fhngz = 0
    end
    fhngs = (90 - del) / 45

    fhng = c10 * fhngrx * fhngr * fhngm * fhngz * fhngs
    # basin response term
    if (z25 <= 1)
      fsed = (c14 + c15 * sj) * (z25 - 1)
    elsif (z25 > 1 and z25 <= 3)
      fsed = 0
    else
      fsed = c16 * k3 * Math.exp(-0.75) * (1 - Math.exp(-0.25 * (z25 - 3)))
    end
    # hypocentral depth term
    if (zhyp <= 7)
      fhyph = 0
    elsif (zhyp > 7 and zhyp <= 20)
      fhyph = zhyp - 7
    else
      fhyph = 13
    end
    if (m <= 5.5)
      fhypm = c17
    elsif (m > 5.5 and m <= 6.5)
      fhypm = c17 + (c18 - c17) * (m - 5.5)
    else
      fhypm = c18
    end
    fhyp = fhyph * fhypm
    # rupture dip term
    if (m <= 4.5)
      fdip = c19 * del
    elsif (m > 4.5 and m <= 5.5)
      fdip = c19 * (5.5 - m) * del
    else
      fdip = 0
    end
    # anelastic attenuation term
    if (region == 0) # ca
      dc20 = dc20ca
    elsif (region == 1 or region == 4) # jp and italy
      dc20 = dc20jp
    elsif (region == 3) # ch
      dc20 = dc20ch
    else
      dc20 = dc20ca # global
    end
    case
      when rrup > 80
        fatn = (c20 + dc20) * (rrup - 80)
      else
        fatn = 0
    end
    if (z25r < 1)
      a1100 = c0 + fmag + fdis + fflt + fhng + fhyp + fdip + fatn + (c11 + k2 * n) * Math.log(1100 / k1) + sj * (c13 + k2 * n) * Math.log(1100 / k1) + (c14 + sj * c15) * (z25r - 1)
    elsif (z25r > 3)
      a1100 = c0 + fmag + fdis + fflt + fhng + fhyp + fdip + fatn + (c11 + k2 * n) * Math.log(1100 / k1) + sj * (c13 + k2 * n) * Math.log(1100 / k1) + c16 * k3 * Math.exp(-0.75) * (1 - Math.exp(-0.25 * (z25r - 3)))
    else
      a1100 = c0 + fmag + fdis + fflt + fhng + fhyp + fdip + fatn + (c11 + k2 * n) * Math.log(1100 / k1) + sj * (c13 + k2 * n) * Math.log(1100 / k1)
    end
    ## dim y as double
    y = Math.exp(a1100)
    a1100_cb = y

    return a1100_cb
  end




  def gmpe_cb_14(m,rrup,rjb,rx,frv,fnm,fhw,ztor,w,del,vs30,z25,zhyp,ztord,wd,zhypd,region,a,t,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,dc20ca,dc20jp,dc20ch,a2,h1,h2,h3,h4,h5,h6,k1,k2,k3,c,n,phi1,phi2,tau1,tau2,philnaf,phic,ro,ypgacb)

    # ***************************************************************
    # *                                                             *
    # *        campbell-bozorgnia nga model, eqs version            *
    # *               written by emel seyhan, phd                   *
    # *                  last checked 03 23 14                      *
    # *                                                             *
    # ***************************************************************
    # this populates the arrays so we can use the correct coefficients when
    #    calculating the formula
    # the period is derived from t and the constants start from
    #    c0 being cof(0) to rlnpga,lny being cof(44)
    ## dim counter as integer
    # counter = 0
    # ## dim cof(43) as double
    # case t
    pi = 4 * Math.atan(1)
    if (region == 1) # jp
      sj = 1
    else
      sj = 0
    end
    # default ztor
    if (ztor == 999 || ztor == nil)
      ztor = ztord
    else
      ztor = ztor
    end
    # default hypocentral depth
    if (zhyp == 999)
      zhyp = zhypd
    else
      zhyp = zhyp
    end
    if (w == 999)
      zhyp = 9
    else
      zhyp = zhyp
    end
    # default w
    if (w == 999)
      w = wd
    else
      w = w
    end
    # magnitude scaling term
    case
      when m <= 4.5
        fmag = (c1 * m)
      when m<= 5.5
        fmag = (c1 * m) + (c2 * (m - 4.5))
      when m<= 6.5
        fmag = (c1 * m) + (c2 * (m - 4.5)) + (c3 * (m - 5.5))
      else
        fmag = (c1 * m) + (c2 * (m - 4.5)) + (c3 * (m - 5.5)) + (c4 * (m - 6.5))
    end
    # distance scaling term
    fdis = ((c5 + (c6 * m)) * (Math.log(Math.sqrt((rrup * rrup) + (c7 * c7)))))
    # style of faulting term
    ffltf = (c8 * frv) + (c9 * fnm)
    case
      when m<= 4.5
        ffltm = 0
      when m<= 5.5
        ffltm = m - 4.5
      else
        ffltm = 1
    end
    fflt = ffltf * ffltm
    # hanging wall term
    r1 = w * Math.cos(del * pi / 180)
    r2 = 62 * m - 350
    f1_rx = h1 + h2 * (rx / r1) + h3 * (rx / r1)**2
    f2_rx = h4 + h5 * ((rx - r1) / (r2 - r1)) + h6 * ((rx - r1) / (r2 - r1))**2
    if (f2_rx > 0)
      maxr = f2_rx
    else
      maxr = 0
    end
    if (fhw == 0)
      fhngrx = 0
    elsif (rx < r1 and fhw == 1)
      fhngrx = f1_rx
    elsif (rx >= r1 and fhw == 1)
      fhngrx = maxr
    end
    case
      when rrup= 0
        fhngr = 1
      else
        fhngr = (rrup - rjb) / rrup
    end
    if (f2_rx > 0)
      f2x = f2_rx
    else
      f2x = 0
    end
    case
      when m<= 5.5
        fhngm = 0
      when m<= 6.5
        fhngm = (m - 5.5) * (1 + a2 * (m - 6.5))
      else
        fhngm = 1 + a2 * (m - 6.5)
    end
    case
      when ztor< 16.66
        fhngz = 1 - 0.06 * ztor
      else
        fhngz = 0
    end
    fhngs = (90 - del) / 45
    fhng = c10 * fhngrx * fhngr * fhngm * fhngz * fhngs
    # shallow site response term
    if (region == 1) # jp
      sjp = 1
    else
      sjp = 0 # others
    end
    if (vs30 <= k1)
      fsite = c11 * Math.log(vs30 / k1) + k2 * (Math.log(a + c * (vs30 / k1)**n) - Math.log(a + c))
    else
      fsite = (c11 + k2 * n) * Math.log(vs30 / k1)
    end
    if (vs30 <= 200)
      fsitej = (c12 + k2 * n) * (Math.log(vs30 / k1) - Math.log(200 / k1)) * sjp
    else
      fsitej = (c13 + k2 * n) * Math.log(vs30 / k1) * sjp
    end
    fsite = fsite + fsitej
    if (z25 == 999)
      if (region == 1) # jp
        z25 = Math.exp(5.359 - 1.102 * Math.log(1100))
      else
        z25 = Math.exp(7.089 - 1.144 * Math.log(1100)) # ca, others
      end
    end
    # basin response term
    if (z25 <= 1)
      fsed = (c14 + c15 * sj) * (z25 - 1)
    elsif (z25 > 1 and z25 <= 3)
      fsed = 0
    else
      fsed = c16 * k3 * Math.exp(-0.75) * (1 - Math.exp(-0.25 * (z25 - 3)))
    end
    # hypocentral depth term
    if (zhyp <= 7)
      fhyph = 0
    elsif (zhyp > 7 and zhyp <= 20)
      fhyph = zhyp - 7
    else
      fhyph = 13
    end
    if (m <= 5.5)
      fhypm = c17
    elsif (m > 5.5 and m <= 6.5)
      fhypm = c17 + (c18 - c17) * (m - 5.5)
    else
      fhypm = c18
    end
    fhyp = fhyph * fhypm
    # rupture dip term
    if (m <= 4.5)
      fdip = c19 * del
    elsif (m > 4.5 and m <= 5.5)
      fdip = c19 * (5.5 - m) * del
    else
      fdip = 0
    end
    # anelastic attenuation term
    if (region == 0) # global, ca, nz, tw
      dc20 = dc20ca
    elsif (region == 1 or region == 4) # jp and italy
      dc20 = dc20jp
    elsif (region == 3) # ch
      dc20 = dc20ch
    else
      dc20 = dc20ca
    end
    case
      when rrup> 80
        fatn = (c20 + dc20) * (rrup - 80)
      else
        fatn = 0
    end
    # model prediction in ln units
    ## dim y as double
    median = Math.exp(c0 + fmag + fdis + fflt + fhng + fsite + fsed + fhyp + fdip + fatn)
    if (t == 0)
      ypgacb = median
    end
    if (t < 0.25 and median < ypgacb)
      y = ypgacb
    else
      y = median
    end
    cb_14 = y

    return cb_14
  end


  def gmpe_cb_14stdev(m,rrup,rjb,rx,frv,fnm,fhw,ztor,w,del,vs30,z25,zhyp,ztord,wd,zhypd,region,a1100,t,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,dc20ca,dc20jp,dc20ch,a2,h1,h2,h3,h4,h5,h6,k1,k2,k3,c,n,phi1,phi2,tau1,tau2,philnaf,phic,ro)

    # ***************************************************************
    # *                                                             *
    # *        campbell-bozorgnia nga model, eqs version            *
    # *               written by silvia mazzoni, phd                   *
    # *                  last checked 03 24 14                      *
    # *                                                             *
    # ***************************************************************

    philnaf0 = 0.3
    phi10 = 0.734
    phi20 = 0.429
    tau10 = 0.409
    tau20 = 0.322

    if (vs30<k1)
      alpha =k2*a1100*((a1100+c*(vs30/k1)**n)**(-1)-(a1100+c)**(-1))
    else
      alpha =0
    end
    if (m<4.5)
      philogy = phi1
    elsif (m>5.5)
      philogy = phi2
    else
      philogy = phi2+(phi1-phi2)*(5.5-m)*(5.5-4.5)
    end
    if (m<4.5)
      philogy0 = phi10
    elsif (m>5.5)
      philogy0 = phi20
    else
      philogy0 = phi20+(phi10-phi20)*(5.5-m)*(5.5-4.5)
    end

    phi=Math.sqrt((philogy**2-philnaf**2)+philnaf**2+alpha**2*(philogy0**2-philnaf0 **2)+2*alpha*ro*Math.sqrt(philogy**2-philnaf**2)*Math.sqrt(philogy0**2-philnaf0 **2))

    if (m<4.5)
      taulogy = tau1
    elsif (m>5.5)
      taulogy = tau2
    else
      taulogy = tau2+(tau1-tau2)*(5.5-m)*(5.5-4.5)
    end

    if (m<4.5)
      taulogy0 = tau10
    elsif(m>5.5)
      taulogy0 = tau20
    else
      taulogy0 = tau20+(tau10-tau20)*(5.5-m)*(5.5-4.5)
    end
    tau = Math.sqrt(taulogy**2+alpha**2*taulogy0**2+2*alpha*ro*taulogy*taulogy0)
    sigma = Math.sqrt(phi**2+tau**2)
    sigmaarb = Math.sqrt(phic**2+sigma**2)
    stdev_cb14 = sigma

  end
  ## public ypga as double
  def gmpe_cy_14(m,rrup,rjb,rx,vs30,frv,fnm,fhw,del,ztor,region,z1,z1r,ddpp,t,c2,c4,c4a,crb,c8,c8a,c1,c1a,c1b,c1c,c1d,cn,cm,c3,c5,chm,c6,c7,c7b,c8b,c9,c9a,c9b,c11,c11b,cg1,cg2,cg3,f1,f2,f3,f4,f5,f6,t1,t2,s1,s2,s3,s2jp,gjpit,gwn,f1jp,f5jp,f6jp,ypga)

    # ***************************************************************
    # *                                                             *
    # *          chiou-youngs nga model, eqs version                *
    # *               written by emel seyhan, phd                   *
    # *                  last checked 03 23 14                      *
    # *                                                             *
    # ***************************************************************
    # this populates the arrays so we can use the correct coefficients when
    #    calculating the formula
    # the period is derived from t and the constants start from
    #    c0 being cof(0) to rlnpga,lny being cof(45)
    ## dim counter as integer
    # counter = 0
    # ## dim cof(45) as double
    # case t
    # when = 0.01
    #     for each cell in range("cy14_coeffs!b5:cy14_coeffs!at5")
    #         cof(counter) = cell.value
    #         counter = counter + 1
    #     next cell
    #  when = 0.02
    #  for each cell in range("cy14_coeffs!b6:cy14_coeffs!at6")
    #         cof(counter) = cell.value
    #         counter = counter + 1
    #     next cell
    #  when = 0.03
    # gjpit = cof(40)
    # gwn = cof(41)
    # f1jp = cof(42)
    # f5jp = cof(43)
    # f6jp = cof(44)
    pi = 4 * Math.atan(1)
    # region: 0 (global, ca, taiwan, new zealand), 1 (jp), 3 (ch), 4 (italy), 5 (tur)
    if (region == 1)
      f1 = f1jp
    end
    if (region == 1)
      s2 = s2jp
    end
    if (region == 1)
      f5 = f5jp
    end
    if (region == 1 or region == 4) # jp and italy
      g = gjpit
    end
    if (region == 3) # ch
      g = gwn
    end
    if (region == 1) # jp
      f6 = f6jp
    end
    # magnitude scaling term
    fmag = c2 * (m - 6) + (c2 - c3) / cn * Math.log(1 + Math.exp(cn * (cm - m)))
    # distance scaling term
    if ((m - chm) > 0)
      maxoftwo = (m - chm)
    else
      maxoftwo = 0
    end
    x = c6 * maxoftwo
    # distance attenuation term
    if ((m - cg3) > 0)
      e = m - cg3
    else
      e = 0
    end
    if ((region == 1 or region == 4) and (6 < m and m < 6.9)) # jp and italy
      fatn = gjpit * (cg1 + cg2 / ((Math.exp(e) + Math.exp(-e)) / 2)) * rrup
    elsif (region == 3)
      fatn = gwn * (cg1 + cg2 / ((Math.exp(e) + Math.exp(-e)) / 2)) * rrup
    else
      fatn = (cg1 + cg2 / ((Math.exp(e) + Math.exp(-e)) / 2)) * rrup
    end
    fdis = c4 * Math.log(rrup + c5 * (Math.exp(x) + Math.exp(-x)) / 2) + (c4a - c4) * Math.log(Math.sqrt(rrup**2 + crb**2)) + fatn
    # style of faulting term
    if ((m - 4.5) > 0)
      maxoftwo2 = m - 4.5
    else
      maxoftwo2 = 0
    end
    z = 2 * maxoftwo2
    fflt = frv * (c1a + c1c / ((Math.exp(z) + Math.exp(-z)) / 2)) + fnm * (c1b + c1d / ((Math.exp(z) + Math.exp(-z)) / 2))
    # hanging wall term
    k = rx / c9b
    tanh = (Math.exp(k) - Math.exp(-k)) / (Math.exp(k) + Math.exp(-k))
    # ztor term
    # center z_tor on the z_tor-m relation in chiou and youngs (2013)
    if (2.704 - 1.226 * (m - 5.849) > 0) # works
      maxztorm = 2.704 - 1.226 * (m - 5.849)
    else
      maxztorm = 0
    end
    if (m <= 5.849 and frv == 1) # works
      mz_torx = 2.704 * 2.704
    else
      mz_torx = maxztorm * maxztorm
    end
    if (2.673 - 1.136 * (m - 4.97) > 0 and fnm == 1)
      maxztorz = 2.673 - 1.136 * (m - 4.97)
    else
      maxztorz = 0
    end
    if (m <= 4.97 and fnm == 1)
      mz_torz = 2.673 * 2.673
    else
      mz_torz = maxztorz * maxztorz
    end
    if (2.673 - 1.136 * (m - 4.97) > 0 and fnm == 0 and frv == 0)
      maxztord = 2.673 - 1.136 * (m - 4.97)
    else
      maxztord = 0
    end
    if (frv == 1)
      mz_tor = mz_torx
    elsif (fnm == 1)
      mz_tor = mz_torz
    else
      mz_tor = maxztord * maxztord
    end
    if (ztor == 999) # default: ca model
      ztor = mz_tor
    else
      ztor = ztor
    end
    deltaz_tor = ztor - mz_tor
    ftor = (c7 + c7b / ((Math.exp(z) + Math.exp(-z)) / 2)) * deltaz_tor
    if (fhw == 0)
      fhng = 0
    else
      fhng = c9 * fhw * Math.cos(del * pi / 180) * (c9a + (1 - c9a) * tanh) * (1 - Math.sqrt(rjb**2 + ztor**2) / (rrup + 1))
    end
    # basin depth term
    # center z1 on the z1-m relation in chiou and youngs (2013)
    if (z1 == 999)
      z1 = z1r
    else
      z1 = z1
    end
    if (region == 1) # jp model
      mz1 = z1 * 1000 - Math.exp(-5.23 / 2 * Math.log((vs30**2 + 412.39**2) / (1360**2 + 412.39**2)))
    else # use default ca model for non-japanese regions
      mz1 = z1 * 1000 - Math.exp(-7.15 / 4 * Math.log((vs30**4 + 570.94**4) / (1360**4 + 570.94**4)))
    end
    if (z1 = z1r)
      deltaz1 = 0
    else
      deltaz1 = mz1
    end
    # dip term
    fdip = (Math.cos(del * pi / 180)**2) * (c11 + c11b / ((Math.exp(z) + Math.exp(-z)) / 2))
    # directivity term
    if ((rrup - 40) > 0)
      maxdir = rrup - 40
    else
      maxdir = 0
    end
    if (1 - maxdir / 30 > 0)
      maxdir1 = 1 - maxdir / 30
    else
      maxdir1 = 0
    end
    if ((m - 5.5) > 0)
      maxdir2 = m - 5.5
    else
      maxdir2 = 0
    end
    if (maxdir2 / 0.8 < 1)
      mindir = maxdir2 / 0.8
    else
      mindir = 1
    end
    fdir = c8 * maxdir1 * mindir * Math.exp(-c8a * (m - c8b)**2) * ddpp
    # yref term
    yref = Math.exp(c1 + fmag + fdis + fflt + fhng + ftor + fdip + fdir)
    # site response term
    if (Math.log(vs30 / 1130) < 0)
      minoftwo = Math.log(vs30 / 1130)
    else
      minoftwo = 0
    end
    if (vs30 < 1130)
      minoftwo2 = vs30
    else
      minoftwo2 = 1130
    end
    rkdepth = f5 * (1 - Math.exp(-deltaz1 / f6))
    fsite = f1 * minoftwo + f2 * (Math.exp(f3 * (minoftwo2 - 360)) - Math.exp(f3 * (1130 - 360))) * Math.log((yref + f4) / f4)
    # model prediction in ln units
    ## dim y as double
    # y term
    y = yref * Math.exp(fsite + rkdepth)
    if (t == 0)
      ypga = y
    end
    if (t <= 0.3 and y < ypga)
      y = ypga
    else
      y = y
    end
    cy_14 = y

    return cy_14
  end




  # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def gmpe_cy14_stdev(m,rrup,rjb,rx,vs30,frv,fnm,fhw,del,ztor,region,z1,ddpp,vs30flag,t,c2,c4,c4a,crb,c8,c8a,c1,c1a,c1b,c1c,c1d,cn,cm,c3,c5,chm,c6,c7,c7b,c8b,c9,c9a,c9b,c11,c11b,cg1,cg2,cg3,f1,f2,f3,f4,f5,f6,t1,t2,s1,s2,s3,s2jp,gjpit,gwn,f1jp,f5jp,f6jp)

    # ***************************************************************
    # *                                                             *
    # *             chiou-youngs sigma, eqs version                 *
    # *               written by emel seyhan, phd                   *
    # *                  last checked 03 23 14                      *
    # *                                                             *
    # ***************************************************************
    ## dim counter as integer
    # counter = 0
    # ## dim cof(45) as double
    # case t
    # when = 0.01
    #     for each cell in range("cy14_coeffs!b5:cy14_coeffs!at5")
    #         cof(counter) = cell.value
    #         counter = counter + 1
    #     next cell
    #  when = 0.02
    # s2jp = cof(39)
    # gjpit = cof(40)
    # gwn = cof(41)
    # f1jp = cof(42)
    # f5jp = cof(43)
    # f6jp = cof(44)
    pi = 4 * Math.atan(1)
    # region: 0 (global, ca, taiwan, new zealand), 1 (jp), 3 (ch), 4 (italy), 5 (tur)
    if (region == 1) # jp
      f1 = f1jp
    end
    if (region == 1) # jp
      s2 = s2jp
    end
    if (region == 1) # jp
      f5 = f5jp
    end
    if (region == 1 or region == 4) # jp and italy
      g = gjpit
    end
    if (region == 3) # ch
      g = gwn
    end
    if (region == 1) # jp
      f6 = f6jp
    end
    # magnitude scaling term
    fmag = c2 * (m - 6) + (c2 - c3) / cn * Math.log(1 + Math.exp(cn * (cm - m)))
    # distance scaling term
    if ((m - chm) > 0)
      maxoftwo = (m - chm)
    else
      maxoftwo = 0
    end
    x = c6 * maxoftwo
    # distance attenuation term
    if ((m - cg3) > 0)
      e = m - cg3
    else
      e = 0
    end
    if ((region == 1 or region == 4) and (6 < m and m < 6.9)) # jp and italy
      fatn = gjpit * (cg1 + cg2 / ((Math.exp(e) + Math.exp(-e)) / 2)) * rrup
    elsif (region == 3)
      fatn = gwn * (cg1 + cg2 / ((Math.exp(e) + Math.exp(-e)) / 2)) * rrup
    else
      fatn = (cg1 + cg2 / ((Math.exp(e) + Math.exp(-e)) / 2)) * rrup
    end
    fdis = c4 * Math.log(rrup + c5 * (Math.exp(x) + Math.exp(-x)) / 2) + (c4a - c4) * Math.log(Math.sqrt(rrup**2 + crb**2)) + fatn
    # style of faulting term
    if ((m - 4.5) > 0)
      maxoftwo2 = m - 4.5
    else
      maxoftwo2 = 0
    end
    z = 2 * maxoftwo2
    fflt = frv * (c1a + c1c / ((Math.exp(z) + Math.exp(-z)) / 2)) + fnm * (c1b + c1d / ((Math.exp(z) + Math.exp(-z)) / 2))
    # hanging wall term
    k = rx / c9b
    tanh = (Math.exp(k) - Math.exp(-k)) / (Math.exp(k) + Math.exp(-k))
    # ztor term
    # center z_tor on the z_tor-m relation in chiou and youngs (2013)
    if (2.704 - 1.226 * (m - 5.849) > 0) # works
      maxztorm = 2.704 - 1.226 * (m - 5.849)
    else
      maxztorm = 0
    end
    if (m <= 5.849 and frv == 1) # works
      mz_torx = 2.704 * 2.704
    else
      mz_torx = maxztorm * maxztorm
    end
    if (2.673 - 1.136 * (m - 4.97) > 0 and fnm == 1)
      maxztorz = 2.673 - 1.136 * (m - 4.97)
    else
      maxztorz = 0
    end
    if (m <= 4.97 and fnm == 1)
      mz_torz = 2.673 * 2.673
    else
      mz_torz = maxztorz * maxztorz
    end
    if (2.673 - 1.136 * (m - 4.97) > 0 and fnm == 0 and frv == 0)
      maxztord = 2.673 - 1.136 * (m - 4.97)
    else
      maxztord = 0
    end
    if (frv == 1)
      mz_tor = mz_torx
    elsif (fnm == 1)
      mz_tor = mz_torz
    else
      mz_tor = maxztord * maxztord
    end
    if (ztor == 999) # default: ca model
      ztor = mz_tor
    else
      ztor = ztor
    end
    deltaz_tor = ztor - mz_tor
    ftor = (c7 + c7b / ((Math.exp(z) + Math.exp(-z)) / 2)) * deltaz_tor
    if (fhw == 0)
      fhng = 0
    else
      fhng = c9 * fhw * Math.cos(del * pi / 180) * (c9a + (1 - c9a) * tanh) * (1 - Math.sqrt(rjb**2 + ztor**2) / (rrup + 1))
    end
    # dip term
    fdip = (Math.cos(del * pi / 180)**2) * (c11 + c11b / ((Math.exp(z) + Math.exp(-z)) / 2))
    # directivity term
    if ((rrup - 40) > 0)
      maxdir = rrup - 40
    else
      maxdir = 0
    end
    if (1 - maxdir / 30 > 0)
      maxdir1 = 1 - maxdir / 30
    else
      maxdir1 = 0
    end
    if ((m - 5.5) > 0)
      maxdir2 = m - 5.5
    else
      maxdir2 = 0
    end
    if (maxdir2 / 0.8 < 1)
      mindir = maxdir2 / 0.8
    else
      mindir = 1
    end
    fdir = c8 * maxdir1 * mindir * Math.exp(-c8a * (m - c8b)**2) * ddpp
    # yref term
    yref = Math.exp(c1 + fmag + fdis + fflt + fhng + ftor + fdip + fdir)
    # site response term
    if (Math.log(vs30 / 1130) < 0)
      minoftwo = Math.log(vs30 / 1130)
    else
      minoftwo = 0
    end
    if (vs30 < 1130)
      minoftwo2 = vs30
    else
      minoftwo2 = 1130
    end
    b = f2 * (Math.exp(f3 * (minoftwo2 - 360)) - Math.exp(f3 * (1130 - 360)))
    nl0 = b * (yref / (yref + f4))
    if (m >= 5)
      maxt = m
    else
      maxt = 5
    end
    if (maxt < 6.5)
      mint = maxt
    else
      mint = 6.5
    end
    tau = t1 + (t2 - t1) / 1.5 * (mint - 5)
    sigma_nl0_s = s1 + (s2 - s1) / 1.5 * (mint - 5)
    if (vs30flag == 1)
      fmeas = 1
      finf = 0
    elsif (vs30flag == 0)
      finf = 1
      fmeas = 0
    end
    sigma_nl0 = sigma_nl0_s * Math.sqrt(0.7 * fmeas + finf * s3 + (1 + nl0)**2)
    sigma = Math.sqrt((tau * (1 + nl0))**2 + sigma_nl0**2)
    phi = sigma_nl0
    tau = (tau * (1 + nl0))
    # model prediction in ln units
    ## dim s as double
    # s term
    s = sigma
    cy14_stdev = s

    return cy14_stdev
  end









  # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  #  use this function for unknown z1
  def gmpe_ask14_z1(vs30,region,z1)

    # ***************************************************************
    # *                                                             *
    # *        abrahamson-silva-kamai basin depth, eqs version      *
    # *               written by emel seyhan, phd                   *
    # *                  last checked 03 23 14                      *
    # *                                                             *
    # ***************************************************************
    if (z1 == 999)
      if (region == 1)
        z1 = Math.exp(-5.23 / 2 * Math.log((vs30**2 + 412**2) / (1360**2 + 412**2))) / 1000
      else # non-japanese
        z1 = Math.exp((-7.67 / 4) * Math.log((vs30**4 + 610**4) / (1360**4 + 610**4))) / 1000
      end
    else
      z1 = z1
    end
    ## dim y as double
    y = z1
    ask14_z1 = y

    return ask14_z1
  end




  def gmpe_ask_14(m,rrup,rjb,rx,frv,fnm,fhw,fas,ztor,w,dip,vs30,vs30flag,z1,ry0,region,t,vlin,b,n,m1,c,c4,a1,a2,a3,a4,a5,a6,a7,a8,a10,a11,a12,a13,a14,a15,a16,a17,a43,a44,a45,a46,a25,a28,a29,a31,a36,a37,a38,a39,a40,a41,a42,s1,s2,s3,s4,s1v,s2v,s5,s6)
    # ***************************************************************
    # *                                                             *
    # *        abrahamson-silva-kamai nga model, eqs version        *
    # *               written by emel seyhan, phd                   *
    # *                  last checked 03 23 14                      *
    # *                                                             *
    # ***************************************************************
    # this populates the arrays so we can use the correct coefficients when
    #    calculating the formula
    # the period is derived from t and the constants start from
    #    c0 being cof(0) to rlnpga,lny being cof(45)
    ## dim counter as integer
    # counter = 0
    # ## dim cof(45) as double
    # case t
    # when = 0.01
    #     for each cell in range("ask14_coeffs!b5:ask14_coeffs!at5")
    #         cof(counter) = cell.value
    #         counter = counter + 1
    #     next cell
    # s2v = cof(42)
    # s5 = cof(43)
    # s6 = cof(44)
    a2hw = 0.2
    h1 = 0.25
    h2 = 1.5
    h3 = -0.75
    m2 = 5
    n = 1.5
    a7 = 0
    crjb = 999.9
    pi = 4 * Math.atan(1)
    # depth-to-top of rupture model
    if (ztor == 999)
      ## dim maxm as double, maxz as double
      if (frv == 1)
        # reverse
        maxm = m - 5.849
        if (0  > maxm)
          maxm = 0
          maxz = 2.704 - 1.226 * maxm
        end
        if (0  > maxz)
          maxz = 0
        end
      else
        # strike slip, normal
        maxm = m - 4.97
        if (0  > maxm)
          maxm = 0
          maxz = 2.673 - 1.136 * maxm
        end
        if (0  > maxz)
          maxz = 0
        end
      end
      ztor = maxz**2
    end
    # estimations for w and ztor
    if (18 / Math.sin(dip * pi / 180) < 10**(-1.75 + 0.45 * m))
      west = 18 / Math.sin(dip * pi / 180)
    else
      west = 10**(-1.75 + 0.45 * m)
    end
    if (w == 999)
      w = west
    else
      w = w
    end
    m1z = 5
    m2z = 7.2
    maxval = 7.8
    if (m <= m1z)
      ztord = maxval
    elsif (m <= m2z)
      ztord = maxval - (maxval / (m2 - m1)) * (m - m1)
    else
      ztord = 0
    end
    if (ztor == 999)
      ztor = ztord
    else
      ztor = ztor
    end
    # magnitude and distance scaling terms
    if (m >= 5)
      c4m = c4
    elsif (m >= 4)
      c4m = c4 - (c4 - 1) * (5 - m)
    else
      c4m = 1
    end
    r = Math.sqrt(rrup**2 + c4m**2)
    if (m < m2)
      f1 = a1 + a4 * (m2 - m1) + a8 * (8.5 - m2)**2 + a6 * (m - m2) + a7 * (m - m2)**2 + (a2 + a3 * (m2 - m1)) * Math.log(r) + a17 * rrup
    elsif (m < m1)
      f1 = a1 + a4 * (m - m1) + a8 * (8.5 - m)**2 + (a2 + a3 * (m - m1)) * Math.log(r) + a17 * rrup
    else
      f1 = a1 + a5 * (m - m1) + a8 * (8.5 - m)**2 + (a2 + a3 * (m - m1)) * Math.log(r) + a17 * rrup
    end
    # vs30star term
    if (t >= 3)
      v1 = 800
    elsif (t > 0.5)
      v1 = Math.exp(-0.35 * Math.log(t / 0.5) + Math.log(1500))
    else
      v1 = 1500
    end
    if (vs30 < v1)
      vs30star = vs30
    else
      vs30star = v1
    end
    # depth to top of rupture term
    if (ztor <= 20)
      f6 = a15 * ztor / 20
    elsif (ztor > 20)
      f6 = a15
    end
    # style of faulting term
    if (m < 4)
      f7 = 0
    elsif (m <= 5)
      f7 = a11 * (m - 4)
    else
      f7 = a11
    end
    if (m < 4)
      f8 = 0
    elsif (m <= 5)
      f8 = a12 * (m - 4)
    else
      f8 = a12
    end

    # hanging wall term
    r1 = w * Math.cos(dip * pi / 180)
    r2 = 3 * r1
    ry1 = rx * Math.tan(20 * pi / 180)
    if (dip > 30)
      t1 = (90 - dip) / 45
    else
      t1 = 60 / 45
    end
    if (m > 6.5)
      t2 = 1 + a2hw * (m - 6.5)
    elsif (m > 5.5)
      t2 = 1 + a2hw * (m - 6.5) - (1 - a2hw) * (m - 6.5)**2
    else
      t2 = 0
    end
    if (rx <= r1)
      t3 = h1 + h2 * (rx / r1) + h3 * (rx / r1)**2
    elsif (rx < r2)
      t3 = 1 - (rx - r1) / (r2 - r1)
    else
      t3 = 0
    end
    if (ztor < 10)
      t4 = 1 - ztor**2 / 100
    else
      t4 = 0
    end
    # compute hw taper 5 (e1q.13) ry0 version
    # taper off edge
    if (ry0 != 999)
      if (ry0 < ry1)
        t5 = 1
      elsif ((ry0 - ry1) < 5)
        t5 = 1 - (ry0 - ry1) / 5
      else
        t5 = 0
      end

    else

      # taper no ry0 version
      if (rjb == 0)
        t5 = 1
      elsif (rjb < 30)
        t5 = 1 - rjb / 30
      else
        t5 = 0
      end
    end
    if (fhw == 1)
      f4 = a13 * t1 * t2 * t3 * t4 * t5
    else
      f4 = 0
    end
    # regionalization term
    if (region == 2) #  for tw
      f12 = a31 * Math.log(vs30star / vlin)
    else
      f12 = 0
    end
    # calculation of sa1180 and site term
    if (region == 1) # jp
      br = (a41 + (1180 - 850) * (a42 - a41) / (1150 - 850))
    elsif (region == 2) # tw
      br = a31 * Math.log(1180 / vlin)
    else
      br = 0
    end
    # br = bx
    if (region == 1) # jp
      by = a29
    elsif (region == 3) # ch
      by = a28
    elsif (region == 2) # tw
      by = a25
    else
      by = 0
    end
    bs = by * rrup + br # bs = bz
    if (1180 >= v1)
      vs30star1180 = v1
    else
      vs30star1180 = 1180
    end
    bq = (a10 + b * n) * Math.log(vs30star1180 / vlin)
    sa1180 = Math.exp(f1 + f6 + frv * f7 + fnm * f8 + fhw * f4 + bq + bs)
    # site amplification
    if (vs30 < vlin)
      f5 = a10 * Math.log(vs30star / vlin) - b * Math.log(sa1180 + c) + b * Math.log(sa1180 + c * (vs30star / vlin)**n)
    else
      f5 = (a10 + b * n) * Math.log(vs30star / vlin)
    end
    # regional z1
    if (region == 1) # jp
      z1r = Math.exp(-5.23 / 2 * Math.log((vs30**2 + 412**2) / (1360**2 + 412**2))) / 1000
    else # non-japanese
      z1r = Math.exp(-7.67 / 4 * Math.log((vs30**4 + 610**4) / (1360**4 + 610**4))) / 1000
    end
    # soil depth model(eq 2/10/2014.17), updated 08/01/13
    if (vs30 <= 150)
      y1z = a43
      y2z = a43
      x1z = 50
      x2z = 150
    elsif (vs30 <= 250)
      y1z = a43
      y2z = a44
      x1z = 150
      x2z = 250
    elsif (vs30 <= 400)
      y1z = a44
      y2z = a45
      x1z = 250
      x2z = 400
    elsif (vs30 <= 700)
      y1z = a45
      y2z = a46
      x1z = 400
      x2z = 700
    else
      y1z = a46
      y2z = a46
      x1z = 700
      x2z = 1000
    end
    # f10 term goes to zero at 1180 m/s (reference)
    if (vs30 == 1180)
      f10 = 0
    else

      f10 = (y1z + (vs30 - x1z) * (y2z - y1z) / (x2z - x1z)) * Math.log((z1 + 0.01) / (z1r + 0.01))
    end
    # aftershock scaling term (eq 4.21)
    if (fas == 1 and crjb >= 15)
      f11 = 0
    elsif (crjb <= 5)
      f11 = a14
    else
      f11 = a14 * (1 - (crjb - 5) / 10)
    end
    if (fas == 0)
      f11 = 0
    end
    # regionalization term for vs30 scaling
    # region == [2,3,1,0] tw, ch, jp, other
    # for japan
    if (region == 1)
      if (vs30 < 150)
        y1 = a36
        y2 = a36
        x1 = 50
        x2 = 150
      elsif (vs30 < 250)
        y1 = a36
        y2 = a37
        x1 = 150
        x2 = 250
      elsif (vs30 < 350)
        y1 = a37
        y2 = a38
        x1 = 250
        x2 = 350
      elsif (vs30 < 450)
        y1 = a38
        y2 = a39
        x1 = 350
        x2 = 450
      elsif (vs30 < 600)
        y1 = a39
        y2 = a40
        x1 = 450
        x2 = 600
      elsif (vs30 < 850)
        y1 = a40
        y2 = a41
        x1 = 600
        x2 = 850
      elsif (vs30 < 1150)
        y1 = a41
        y2 = a42
        x1 = 850
        x2 = 1150
      else
        y1 = a42
        y2 = a42
        x1 = 1150
        x2 = 3000
      end
      f13 = y1 + (y2 - y1) / (x2 - x1) * (vs30 - x1)
    end
    # ftw = 1 for tw, fcn = 1 for ch, fjp = 1 for jp, 0 for other
    if (region == 2) # tw for ask14
      reg = f12 + a25 * rrup
    elsif (region == 3) # ch
      reg = a28 * rrup
    elsif (region == 1) # jp
      reg = f13 + a29 * rrup
    else
      reg = 0
    end
    # regional vs30 scaling term # bx
    if (region == 1) # jp
      rvs = f13
    elsif (region == 2) # tw
      rvs = f12
    else
      rvs = 0
    end
    # regional delta
    rd = rvs + by * rrup
    # model prediction in ln units
    ## dim y as double

    y = Math.exp(f1 + frv * f7 + fnm * f8 + fas * f11 + f5 + fhw * f4 + f6 + f10 + rd)
    ask_14 = y

    return ask_14
  end




  # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def gmpe_ask14_stdev(m,rrup,rjb,rx,frv,fnm,fhw,fas,ztor,w,dip,vs30,vs30flag,z1,ry0,region,t,vlin,b,n,m1,c,c4,a1,a2,a3,a4,a5,a6,a7,a8,a10,a11,a12,a13,a14,a15,a16,a17,a43,a44,a45,a46,a25,a28,a29,a31,a36,a37,a38,a39,a40,a41,a42,s1,s2,s3,s4,s1v,s2v,s5,s6)
    # ***************************************************************
    # *                                                             *
    # *        abrahamson-silva-kamai sigma, eqs version            *
    # *               written by emel seyhan, phd                   *
    # *                  last checked 03 23 14                      *
    # *                                                             *
    # ***************************************************************
    ## dim counter as integer
    # counter = 0
    # ## dim cof(45) as double
    # case t
    # when = 0.01
    a2hw = 0.2
    h1 = 0.25
    h2 = 1.5
    h3 = -0.75
    m2 = 5
    n = 1.5
    a7 = 0
    samp = 0.4
    crjb = 999.9
    pi = 4 * Math.atan(1)
    # depth-to-top of rupture model
    if (ztor == 999)
      ## dim maxm as double, maxz as double
      if (frv == 1)
        # reverse
        maxm = m - 5.849
        if (0  > maxm)
          maxm = 0
          maxz = 2.704 - 1.226 * maxm
        end
        if (0  > maxz)
          maxz = 0
        end
      else
        # strike slip, normal
        maxm = m - 4.97
        if (0  > maxm)
          maxm = 0
          maxz = 2.673 - 1.136 * maxm
        end
        if (0  > maxz)
          maxz = 0
        end
      end
      ztor = maxz**2
    end
    # estimations for w and ztor
    if (18 / Math.sin(dip * pi / 180) < 10**(-1.75 + 0.45 * m))
      west = 18 / Math.sin(dip * pi / 180)
    else
      west = 10**(-1.75 + 0.45 * m)
    end
    if (w == 999)
      w = west
    else
      w = w
    end
    m1z = 5
    m2z = 7.2
    maxval = 7.8
    if (m <= m1z)
      ztord = maxval
    elsif (m <= m2z)
      ztord = maxval - (maxval / (m2 - m1)) * (m - m1)
    else
      ztord = 0
    end
    if (ztor == 999)
      ztor = ztord
    else
      ztor = ztor
    end
    # magnitude and distance scaling terms
    if (m >= 5)
      c4m = c4
    elsif (m >= 4)
      c4m = c4 - (c4 - 1) * (5 - m)
    else
      c4m = 1
    end
    r = Math.sqrt(rrup**2 + c4m**2)
    if (m < m2)
      f1 = a1 + a4 * (m2 - m1) + a8 * (8.5 - m2)**2 + a6 * (m - m2) + a7 * (m - m2)**2 + (a2 + a3 * (m2 - m1)) * Math.log(r) + a17 * rrup
    elsif (m < m1)
      f1 = a1 + a4 * (m - m1) + a8 * (8.5 - m)**2 + (a2 + a3 * (m - m1)) * Math.log(r) + a17 * rrup
    else
      f1 = a1 + a5 * (m - m1) + a8 * (8.5 - m)**2 + (a2 + a3 * (m - m1)) * Math.log(r) + a17 * rrup
    end
    # vs30star term
    if (t >= 3)
      v1 = 800
    elsif (t > 0.5)
      v1 = Math.exp(-0.35 * Math.log(t / 0.5) + Math.log(1500))
    else
      v1 = 1500
    end
    if (vs30 < v1)
      vs30star = vs30
    else
      vs30star = v1
    end
    # depth to top of rupture term
    if (ztor <= 20)
      f6 = a15 * ztor / 20
    elsif (ztor > 20)
      f6 = a15
    end
    # style of faulting term
    if (m < 4)
      f7 = 0
    elsif (m <= 5)
      f7 = a11 * (m - 4)
    else
      f7 = a11
    end
    if (m < 4)
      f8 = 0
    elsif (m <= 5)
      f8 = a12 * (m - 4)
    else
      f8 = a12
    end
    # hanging wall term
    r1 = w * Math.cos(dip * pi / 180)
    r2 = 3 * r1
    ry1 = rx * Math.tan(20 * pi / 180)
    if (dip > 30)
      t1 = (90 - dip) / 45
    else
      t1 = 60 / 45
    end
    if (m > 6.5)
      t2 = 1 + a2hw * (m - 6.5)
    elsif (m > 5.5)
      t2 = 1 + a2hw * (m - 6.5) - (1 - a2hw) * (m - 6.5)**2
    else
      t2 = 0
    end
    if (rx <= r1)
      t3 = h1 + h2 * (rx / r1) + h3 * (rx / r1)**2
    elsif (rx < r2)
      t3 = 1 - (rx - r1) / (r2 - r1)
    else
      t3 = 0
    end
    if (ztor < 10)
      t4 = 1 - ztor**2 / 100
    else
      t4 = 0
    end
    # compute hw taper 5 (e1q.13) ry0 version
    # taper off edge
    if (ry0 != 999)
      if (ry0 < ry1)
        t5 = 1
      elsif ((ry0 - ry1) < 5)
        t5 = 1 - (ry0 - ry1) / 5
      else
        t5 = 0
      end

    else

      # no ry0 version
      # taper no ry0
      if (rjb == 0)
        t5 = 1
      elsif (rjb < 30)
        t5 = 1 - rjb / 30
      else
        t5 = 0
      end
    end
    if (fhw == 1)
      f4 = a13 * t1 * t2 * t3 * t4 * t5
    else
      f4 = 0
    end
    # calculation of sa1180 and site term
    if (region == 3) # ch
      br = (a41 + (1180 - 850) * (a42 - a41) / (1150 - 850))
    elsif (region == 2) # tw
      br = a31 * Math.log(1180 / vlin)
    else
      br = 0
    end
    # br = bx
    if (region == 1)  # jp
      by = a29
    elsif (region == 3) # ch
      by = a28
    elsif (region == 2) # tw
      by = a25
    else
      by = 0
    end
    bs = by * rrup + br # bs = bz
    if (1180 >= v1)
      vs30star1180 = v1
    else
      vs30star1180 = 1180
    end
    bq = (a10 + b * n) * Math.log(vs30star1180 / vlin)
    sa1180 = Math.exp(f1 + f6 + frv * f7 + fnm * f8 + fhw * f4 + bq + bs)
    if (m < 4)
      phia = s1
    elsif (m < 6)
      phia = s1 + (s2 - s1) * (m - 4) / 2
    else
      phia = s2
    end
    if (region == 1) # jp
      if (rrup < 30)
        phia2 = s5
      elsif (rrup >= 30 and rrup <= 80)
        phia2 = s5 + (s6 - s5) * (rrup - 30) / 50
      elsif (rrup > 80)
        phia2 = s6
      end
    else
      phia2 = phia
    end
    # known vs30
    if (m < 4)
      phiat = s1v
    elsif (m < 6)
      phiat = s1v + (s2v - s1v) * (m - 4) / 2
    else
      phiat = s2v
    end
    if (m < 5)
      taua = s3
    elsif (m < 7)
      taua = s3 + (s4 - s3) * (m - 5) / 2
    else
      taua = s4
    end
    if (vs30star >= vlin)
      ratio = 0
    else
      ratio = b * sa1180 * (-1 / (sa1180 + c) + 1 / (sa1180 + c * (vs30star / vlin)**n))
    end
    tau = taua * (1 + ratio)
    if (vs30flag == 0)
      if (samp**2>phia2**2)
        phib = 0
      else
        phib = Math.sqrt(phia2**2 - samp**2)
      end
    else
      if (samp**2>phiat**2)
        phib = 0
      else
        phib = Math.sqrt(phiat**2 - samp**2)
      end
    end
    phi = Math.sqrt(phib**2 * (1 + ratio)**2 + samp**2)
    sigma = Math.sqrt(tau**2 + phi**2)
    ## dim sd as double
    sd = sigma
    ask14_stdev = sd

    return ask14_stdev
  end




  def gmpe_bssa_14(m,rjb,u,rs,ns,vs30,region,z1,pgar,t,e0,e1,e2,e3,e4,e5,e6,mh,c1,c2,c3,mref,rref,h,dc3,dc3chtur,dc3jpit,c,vc,vref,f1,f3,f4,f5,f6,f7,r1,r2,dfr,dfv,v1,v2,phi1,phi2,tau1,tau2)

    # # *******************************************************************
    # *                                                                  *
    # *    boore-atkinson-stewart-seyhan nga model, version june 2013    *
    # *                   written by emel seyhan, phd                    *
    # *                     last checked 02 10 14                        *
    # *                                                                  *
    # ********************************************************************
    # this populates the arrays so we can use the correct coefficients when
    #    calculating the formula
    # the period is derived from t and the constants start from
    #    c0 being cof(0) to rlnpga,lny being cof(36)
    ## dim counter as integer
    # counter = 0
    # ## dim cof(36) as double
    # case t
    # when
    #     next cell
    #  when = -1
    #     for each cell in range("bssa14_coeffs!b27:bssa14_coeffs!ak27")
    #         cof(counter) = cell.value
    #         counter = counter + 1
    #     next cell
    # end
    # e0 = cof(0)
    # e1 = cof(1)
    # e2 = cof(2)
    # e3 = cof(3)
    # e4 = cof(4)
    # e5 = cof(5)
    # e6 = cof(6)
    # mh = cof(7)
    # c1 = cof(8)
    # c2 = cof(9)
    # c3 = cof(10)
    # mref = cof(11)
    # rref = cof(12)
    # h = cof(13)
    # dc3 = cof(14)
    # dc3chtur = cof(15)
    # dc3jpit = cof(16)
    # c = cof(17)
    # vc = cof(18)
    # vref = cof(19)
    # f1 = cof(20)
    # f3 = cof(21)
    # f4 = cof(22)
    # f5 = cof(23)
    # f6 = cof(24)
    # f7 = cof(25)
    # r1 = cof(26)
    # r2 = cof(27)
    # dfr = cof(28)
    # dfv = cof(29)
    # v1 = cof(30)
    # v2 = cof(31)
    # phi1 = cof(32)
    # phi2 = cof(33)
    # tau1 = cof(34)
    # tau2 = cof(35)
    # magnitude scaling term
    if (ns == 0 and rs == 0 and u == 0)
      ss = 1
    else
      ss = 0
    end
    if (m <= mh)
      fm = e0 * u + e1 * ss + e2 * ns + e3 * rs + e4 * (m - mh) + e5 * (m - mh)**2
    else
      fm = e0 * u + e1 * ss + e2 * ns + e3 * rs + e6 * (m - mh)
    end
    r = Math.sqrt(rjb**2 + h**2)
    if (region == 0) # global, ca, tw, nz
      fp = (c1 + c2 * (m - mref)) * Math.log(r / rref) + (c3 + dc3) * (r - rref)
    elsif (region == 3 or region == 5) # ch, tur
      fp = (c1 + c2 * (m - mref)) * Math.log(r / rref) + (c3 + dc3chtur) * (r - rref)
    elsif (region == 1 or region == 4) # jp, italy
      fp = (c1 + c2 * (m - mref)) * Math.log(r / rref) + (c3 + dc3jpit) * (r - rref)
    else
      fp = (c1 + c2 * (m - mref)) * Math.log(r / rref) + (c3 + dc3) * (r - rref)
    end
    # linear site term
    if (vs30 <= vc)
      flin = c * Math.log(vs30 / vref)
    else
      flin = c * Math.log(vc / vref)
    end
    # nonlinear site term
    if (vs30 < 760)
      minv = vs30
    else
      minv = 760
    end
    f2 = f4 * ((Math.exp(f5 * (minv - 360))) - Math.exp(f5 * (760 - 360)))

    fnl = f1 + f2 * Math.log((pgar + f3) / f3)
    fnl = f1 + (f4 * ((Math.exp(f5 * (minv - 360))) - Math.exp(f5 * (760 - 360)))) * Math.log((pgar + f3) / f3)
    # basin depth term
    if (region == 0) # ca, global, nz, tw
      mz1 = Math.exp(-7.15 / 4 * Math.log((vs30**4 + 570.94**4) / (1360**4 + 570.94**4))) / 1000
    elsif (region == 1) # jp
      mz1 = Math.exp(-5.23 / 2 * Math.log((vs30**2 + 412.39**2) / (1360**2 + 412.39**2))) / 1000
    else # global, default ca
      mz1 = Math.exp(-7.15 / 4 * Math.log((vs30**4 + 570.94**4) / (1360**4 + 570.94**4))) / 1000
    end
    dz1 = z1 - mz1
    if (z1 == 999)
      dz1 = 0
    else
      dz1 = dz1
    end
    if (t < 0.65)
      fz1 = 0
    elsif (dz1 <= f7 / f6)
      fz1 = f6 * dz1
    elsif (dz1 > f7 / f6)
      fz1 = f7
    else
      fz1 = 0
    end
    if (z1 == 999)
      fz1 = 0
    else
      fz1 = fz1
    end
    # site term
    fs = flin + fnl # in ln units
    # model prediction in ln units
    ## dim y as double
    y = Math.exp(fm + fp + fs + fz1)
    bssa_14 = y
    return bssa_14
  end



  # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def gmpe_bssa14_stdev(m,rjb,vs30,t,e0,e1,e2,e3,e4,e5,e6,mh,c1,c2,c3,mref,rref,h,dc3,dc3chtur,dc3jpit,c,vc,vref,f1,f3,f4,f5,f6,f7,r1,r2,dfr,dfv,v1,v2,phi1,phi2,tau1,tau2)

    # # *******************************************************************
    # *                                                                  *
    # *        boore-atkinson-stewart-seyhan sigma, eqs version          *
    # *                   written by emel seyhan, phd                    *
    # *                     last checked 08 12 13                        *
    # *                                                                  *
    # ********************************************************************
    ## dim counter as integer
    # counter = 0
    # ## dim cof(36) as double
    # case t
    # when = 0.01
    #     for each cell in range("bssa14_coeffs!b5:bssa14_coeffs!ak5")
    #         cof(counter) = cell.value
    #         counter = counter + 1
    #     next cell
    #  when = 0.02
    #     next cell
    # end
    # e0 = cof(0)
    # e1 = cof(1)
    # e2 = cof(2)
    # e3 = cof(3)
    # e4 = cof(4)
    # e5 = cof(5)
    # e6 = cof(6)
    # mh = cof(7)
    # c1 = cof(8)
    # c2 = cof(9)
    # c3 = cof(10)
    # mref = cof(11)
    # rref = cof(12)
    # h = cof(13)
    # dc3 = cof(14)
    # dc3chtur = cof(15)
    # dc3jpit = cof(16)
    # c = cof(17)
    # vc = cof(18)
    # vref = cof(19)
    # f1 = cof(20)
    # f3 = cof(21)
    # f4 = cof(22)
    # f5 = cof(23)
    # f6 = cof(24)
    # f7 = cof(25)
    # r1 = cof(26)
    # r2 = cof(27)
    # dfr = cof(28)
    # dfv = cof(29)
    # v1 = cof(30)
    # v2 = cof(31)
    # phi1 = cof(32)
    # phi2 = cof(33)
    # tau1 = cof(34)
    # tau2 = cof(35)
    if (m <= 4.5)
      taum = tau1
    elsif (m > 4.5 and m < 5.5)
      taum = tau1 + (tau2 - tau1) * (m - 4.5)
    else
      taum = tau2
    end
    if (m <= 4.5)
      phim = phi1
    elsif (m > 4.5 and m < 5.5)
      phim = phi1 + (phi2 - phi1) * (m - 4.5)
    else
      phim = phi2
    end
    if (rjb <= r1)
      phimr = phim
    elsif (rjb > r1 and rjb <= r2)
      phimr = phim + dfr * (Math.log(rjb / r1) / (Math.log(r2 / r1)))
    else
      phimr = phim + dfr
    end
    if (vs30 >= v2)
      phi = phimr
    elsif (vs30 >= v1 and vs30 <= v2)
      phi = phimr - dfv * (Math.log(v2 / vs30) / (Math.log(v2 / v1)))
    else
      phi = phimr - dfv
    end
    # model prediction in ln units
    ## dim sigma as double
    sigma = Math.sqrt(taum**2 + phi**2)
    bssa14_stdev = sigma
  end


  # this function is to calculate pgar for bssa14
  def gmpe_pgar_calc(m,rjb,u,rs,ns,region,t,e0,e1,e2,e3,e4,e5,e6,mh,c1,c2,c3,mref,rref,h,dc3,dc3chtur,dc3jpit,c,vc,vref,f1,f3,f4,f5,f6,f7,r1,r2,dfr,dfv,v1,v2,phi1,phi2,tau1,tau2)

    # # *******************************************************************
    # *                                                                  *
    # *        boore-atkinson-stewart-seyhan pgar, eqs version           *
    # *                   written by emel seyhan, phd                    *
    # *                     last checked 03 23 14                        *
    # *                                                                  *
    # ********************************************************************
    ## dim counter as integer
    # counter = 0
    # ## dim cof(17) as double
    # case t
    # when = 0
    #   for each cell in range("bssa14_coeffs!b26:bssa14_coeffs!r26")
    #     cof(counter) = cell.value
    #     counter = counter + 1
    # next cell
    # end
    # e0 = cof(0)
    # e1 = cof(1)
    # e2 = cof(2)
    # e3 = cof(3)
    # e4 = cof(4)
    # e5 = cof(5)
    # e6 = cof(6)
    # mh = cof(7)
    # c1 = cof(8)
    # c2 = cof(9)
    # c3 = cof(10)
    # mref = cof(11)
    # rref = cof(12)
    # h = cof(13)
    # dc3 = cof(14)
    # dc3chtur = cof(15)
    # dc3jpit = cof(16)
    if (ns == 0 and rs == 0 and u == 0)
      ss = 1
    else
      ss = 0
    end
    if (m <= mh)
      fm = e0 * u + e1 * ss + e2 * ns + e3 * rs + e4 * (m - mh) + e5 * (m - mh)**2
    else
      fm = e0 * u + e1 * ss + e2 * ns + e3 * rs + e6 * (m - mh)
    end
    r = Math.sqrt(rjb**2 + h**2)
    if (region == 0) # global, ca, tw, nz
      fp = (c1 + c2 * (m - mref)) * Math.log(r / rref) + (c3 + dc3) * (r - rref)
    elsif (region == 3 or region == 5) # china, turkey
      fp = (c1 + c2 * (m - mref)) * Math.log(r / rref) + (c3 + dc3chtur) * (r - rref)
    elsif (region == 1 or region == 4) # jp, italy
      fp = (c1 + c2 * (m - mref)) * Math.log(r / rref) + (c3 + dc3jpit) * (r - rref)
    else
      fp = (c1 + c2 * (m - mref)) * Math.log(r / rref) + (c3 + dc3) * (r - rref)
    end
    ## dim y as double
    pgar = Math.exp(fm + fp)
    pgar_calc = pgar

    return pgar_calc

  end



  # this code calculates dz1 for bssa14
  def gmpe_dz1_calc(vs30,z1,region)
    # # *******************************************************************
    # *                                                                  *
    # *     boore-atkinson-stewart-seyhan basin depth, eqs version       *
    # *                   written by emel seyhan, phd                    *
    # *                     last checked 03 23 14                        *
    # *                                                                  *
    # ********************************************************************
    # basin depth term
    if (region == 1) # jp
      mz1 = mz1 = Math.exp(-5.23 / 2 * Math.log((vs30**2 + 412.39**2) / (1360**2 + 412.39**2))) / 1000
    else # global
      mz1 = Math.exp(-7.15 / 4 * Math.log((vs30**4 + 570.94**4) / (1360**4 + 570.94**4))) / 1000
    end
    #
    # if (z1 != 999)
    #     dz1 = z1 - mz1
    # elsif (z1 == 999)
    #     dz1 = 0
    # end
    ## dim dz1 as double
    dz1 = z1 - mz1
    dz1_calc = dz1

    return dz1_calc
  end


  def gmpe_i_14(m,rrup,f,vs30,t,a11,a21,b11,b21,a12,a22,b12,b22,a3,ksi,g,phi)
    # # *******************************************************************
    # *                                                                  *
    # *               i.m. idriss nga model, eqs version                 *
    # *                   written by emel seyhan, phd                    *
    # *                     last checked 03 23 14                        *
    # *                                                                  *
    # ********************************************************************
    # this populates the arrays so we can use the correct coefficients when
    #    calculating the formula
    # the period is derived from t and the constants start from
    #    c0 being cof(0) to rlnpga,lny being cof(36)
    ## dim counter as integer
    # counter = 0
    # ## dim cof(12) as double
    # case t
    # when = 0.01
    #         counter = counter + 1
    #     next cell
    # end
    # a11 = cof(0)
    # a21 = cof(1)
    # b11 = cof(2)
    # b21 = cof(3)
    # a12 = cof(4)
    # a22 = cof(5)
    # b12 = cof(6)
    # b22 = cof(7)
    # a3 = cof(8)
    # ksi = cof(9)
    # g = cof(10)
    # phi = cof(11)
    # magnitude scaling term
    if (m <= 6.75)
      fm = a21 * m + a3 * (8.5 - m)**2
    else
      fm = a22 * m + a3 * (8.5 - m)**2
    end
    # distance scaling term
    if (m <= 6.75)
      fdis = -(b11 + b21 * m) * Math.log(rrup + 10)
    else
      fdis = -(b12 + b22 * m) * Math.log(rrup + 10)
    end
    # site term
    if (vs30 < 1200)
      fsite = ksi * Math.log(vs30)
    else
      fsite = ksi * Math.log(1200)
    end
    # attenuation term
    fatn = g * rrup
    # style of faulting term
    fflt = phi * f
    # model prediction in ln units
    ## dim y as double
    if (m <= 6.75)
      y = Math.exp(a11 + fm + fdis + fsite + fatn + fflt)
    else
      y = Math.exp(a12 + fm + fdis + fsite + fatn + fflt)
    end
    i_14 = y

    return i_14
  end


  # # # # # # # # # # # # # # # # # # # # # # # #
  def gmpe_i_14_stdev(m,t,a11,a21,b11,b21,a12,a22,b12,b22,a3,ksi,g,phi)
    # # *******************************************************************
    # *                                                                  *
    # *               i.m. idriss nga model, eqs version                 *
    # *                   written by emel seyhan, phd                    *
    # *                     last checked 03 23 14                        *
    # *                                                                  *
    # ********************************************************************
    # this populates the arrays so we can use the correct coefficients when
    #    calculating the formula
    # the period is derived from t and the constants start from
    #    c0 being cof(0) to rlnpga,lny being cof(36)
    ## dim counter as integer
    # counter = 0
    # ## dim cof(12) as double
    # case t
    # wh
    #     next cell
    # end
    # a11 = cof(0)
    # a21 = cof(1)
    # b11 = cof(2)
    # b21 = cof(3)
    # a12 = cof(4)
    # a22 = cof(5)
    # b12 = cof(6)
    # b22 = cof(7)
    # a3 = cof(8)
    # ksi = cof(9)
    # g = cof(10)
    # phi = cof(11)
    # compute sigma which is period and magnitude dependent.
    # note report does not state a limit on sigma for m<5 but
    # since model is only applicable for m>=5 a limit is retained
    # for sigma with m<5 equal to m=5 values.
    if (t <= 0.05)
      if (m <= 5)
        sigma = 1.18 + 0.035 * Math.log(0.05) - 0.06 * 5
      elsif (m >= 7.5)
        sigma = 1.18 + 0.035 * Math.log(0.05) - 0.06 * 7.5
      else
        sigma = 1.18 + 0.035 * Math.log(0.05) - 0.06 * m
      end
    elsif (t >= 3)
      if (m <= 5)
        sigma = 1.18 + 0.035 * Math.log(3) - 0.06 * 5
      elsif (m >= 7.5)
        sigma = 1.18 + 0.035 * Math.log(3) - 0.06 * 7.5
      else
        sigma = 1.18 + 0.035 * Math.log(3) - 0.06 * m
      end
    else
      if (m <= 5)
        sigma = 1.18 + 0.035 * Math.log(t) - 0.06 * 5
      elsif (m >= 7.5)
        sigma = 1.18 + 0.035 * Math.log(t) - 0.06 * 7.5
      else
        sigma = 1.18 + 0.035 * Math.log(t) - 0.06 * m
      end
    end
    # model prediction
    ## dim y as double
    i_14_stdev = sigma

    return i_14_stdev
  end



  # def AdminFunction_ManageVerifiedUser()
  #   logger.info "in AdminFunction_ManageVerifiedUser"

  #   verifiedEmailList = ["silviamazzoni@yahoo.com","mazzoni@berkeley.edu","silviamazzoni4@yahoo.com","dnewby@absconsulting.com","arahmani@advgeosolutions.com","craig.oldfield@aecom.com","matthew.crake@aecom.com","nik.richter@aecom.com","melanie.walling@aecom.com","patrick.ho2@aecom.com","zavien.teh@aecom.com","mmirmoezi@gmail.com","marianne.malm@afconsult.com","kghiassi@agsinc.com","arnel@aguelassociates.com","mesutturel@gmail.com","aatry@albus-keefe.net","dalbus@albus-keefe.net","bidjan7@hotmail.com","kcallahan@sewardgeo.com","silvia.lohrmann@swissnuclear.ch","hari.ponnaboyina@amec.com","catherine.tatarniuk@amec.com","dixieann.simon@amec.com","alexander.wright@amec.com","eric.mohlmann@amec.com","makram.sabbagh@amec.com","zhihuali1@gmail.com","zhangdyhit@gmail.com","jussi-pekka.matilainen@outlook.com","ashama@ammann-whitney.com","taylorn@agecinc.com","hawkes@agecinc.com","michelle.shriro@agsinc.com","keyvan.fotoohi@agsinc.com","myang@arescorporation.com","erin.leung@arup.com","james.care@arup.com","areti.koskosidi@arup.com","francisco.ciruela-ochoa@arup.com","taki.vlasakakis@arup.com","jongwon.lee@arup.com","malcon.josef@arup.com","rein.de-vries@arup.com","roy.ngan@arup.com","giada.di-fonzo@arup.com","david.dekoning@arup.com","kirk.ellison@arup.com","merrick.taylor@arup.com","jacky.ho@arup.com","gregory.nielsen@arup.com","jorge.lopez@arup.com","andrew.yeskoo@arup.com","minly.so@arup.com","michael.vasileiadis@arup.com","caner.soydas@arup.com","pawan.kumar@arup.com","kendra.jones@arup.com","armin.masroor@arup.com","allan.olson@arup.com","besitler@gmail.com","khalilaa@gmail.com","ayse@atcouncil.org","kedar.kale@atkinsglobal.com","dschwarm@atlasgeotechnical.com","ramin.latifi@ausenco.com","jasonvz@baggengineers.com","brent.bergman@bchydro.com","osmar.penner@bchydro.com","maoxin.li@bchydro.com","sesamory@gmail.com","melanie.regino@beca.com","nick@ngregor.com","huoyili00@gmail.com","eghbal7alireza@gmail.com","rkumar@bhspecialty.com","gomerbm@bv.com","zhengw@bv.com","sommerfeldg@bv.com","jocelyn@bourcetengineering.com","mjbutler@burnsmcd.com","martin.button@sbcglobal.net","michelle.lee@water.ca.gov","mtyler@water.ca.gov","ama@caleng.com","alex@crceng.com","mwoods@water.ca.gov","rarmsto@water.ca.gov","yeo.yoon@dot.ca.gov","seungwoon_han@dot.ca.gov","tom.shantz@dot.ca.gov","cliff.roblee@dot.ca.gov","toorak_zokaie@dot.ca.gov","carlkimgeo@gmail.com","kborromeo@cbi.com","rui.chen@conservation.ca.gov","manuela.davi@ch2m.com","rbethapu@ch2m.com","deepak.rayamajhi@gmail.com","nofalhisham@gmail.com","gabguerr@gmail.com","greg@gemertz.com","epratt@clarkpacific.com","ariyabala@hotmail.com","patxi_uriz@yahoo.com","ccarpenter@cornforthconsultants.com","akost@cornforthconsultants.com","disaacs@crosbygroup.com","alimohammadcsi1990@gmail.com","emilyg@dnfsb.gov","mhachem@degenkolb.com","nacamuli@degenkolb.com","zelzeleh133@yahoo.com","josh.gregorio@dntanks.com","b_raji@hotmail.com","p.wilson@earthmech.com","a.zand@earthmech.com","sstringer9@gmail.com","pespinosa@engeo.com","jfippin@engeo.com","peker@erdemli.com","tmorgan@exponent.com","puriz@exponent.com","adel.h.younan@exxonmobil.com","farzad@fnaeim.com","christopher.wang@ferc.gov","michael.vail@ferc.gov","erin.williams@ferc.gov","vinh.tran@ferc.gov","eric.lim@dot.gov","phayak.takkabutr@fluor.com","jochen.carl@fmglobal.com","m.dana@forell.com","dlr@foundationengr.com","zhiyong.chen@fpinnovations.ca","jve_fredschott@sbcglobal.net","michael.law@fpaengineers.com","chayden@fugro.com","dtsiaousi@fugro.com","umraju@fugro.com","afernandez@fugro.com","glavrentiadis@fugro.com","jkocijan@fugro.com","pwallbridge@fugro.com","wychen@fugro.com","jugalde@fugro.com","d.oconnell@fugro.com","j.altekruse@fugro.com","jmeneses@geiconsultants.com","mh@geocomp.com","thomas@geoconinc.com","arc@geoconsult.us","npaveglio@geodesigninc.com","jhoffman@geoengineers.com","hpuangnak@geoengineers.com","dmclay@geoengineers.com","zsimpson@geoengineers.com","didier.perret@canada.ca","earthquakes@geomotions.com","sismos@geomotions.com","mains@geopacific.ca","doug_wahl@geopentech.com","alexandra_sarmiento@geopentech.com","andrew_dinsick@geopentech.com","carola.dialessandro@gmail.com","jacquelynallmond@gmail.com","zamini@geosyntec.com","dumberg@geosyntec.com","bmartinez@geosyntec.com","grix@geosyntec.com","dadrian25@gmail.com","travis@gerhartcole.com","anirudh.rao@globalquakemodel.org","feng_li@golder.com","cjeong@golder.com","mkennedy@golder.com","aaugello@golder.com","nkoragappa@golder.com","cwoods@golder.com","kbeen@golder.com","mkaya@golder.com","bbayrak@golder.com","nouri_hr@yahoo.com","jbock@gri.com","wspang@gri.com","tmeskele@gri.com","katherineh@groupdelta.com","kristenc@groupdelta.com","bill@gshgeotech.com","kavin@gtcgeotech.com","dustin@gtcgeotech.com","kramerjm@gmail.com","jim.alders@hartcrowser.com","ben.blanchette@hartcrowser.com","chris.delatorre@hartcrowser.com","kevin.huynh@hatchmott.com","marohrbach@haywardbaker.com","bskyers@hdgeoinc.com","josh.zupan@gmail.com","jeffrey.svatora@hdrinc.com","clinton.forsha@hdrinc.com","marriott.dj@gmail.com","srivishnu.mohan@gmail.com","joannagillie@hotmail.com","robert.spears@inl.gov","geotechnicalengineer@gmail.com","ing_anacastillo@hotmail.com","j.lockhart@jmlockhart.net","j.swaisgood@msn.com","greadshaw@jznengineering.com","njundi@jznengineering.com","bdibra1@gmail.com","zzafir@kleinfelder.com","jgingery@kleinfelder.com","yzhou@kleinfelder.com","kulmer@kleinfelder.com","blingwall@kleinfelder.com","s.ainur2010@gmail.com","mpernito@klohn.com","burakkoza@hotmail.com","anna.migliaccio@kpff.com","shane.noel@kpff.com","reid.zimmerman@kpff.com","will.mcvitty@kpff.com","groberts@kpff-la.com","nayche@lafp.com","galcantar@langan.com","jgouchon@langan.com","rgolesorkhi@langan.com","ksyngros@langan.com","tracy.guo.tg@gmail.com","jwilliams@latitudedrc.com","jbwong@lbl.gov","spulijala@leightongroup.com","givler@lettisci.com","jliao@kleinfelder.com","carene@lanl.gov","veronica@maffei-structure.com","r_sojoudi@yahoo.com","hosein.kerdar@marstructuraldesign.com","linzjonz@maxlide.com","taha.ashraf@gmail.com","srashidi@megconsulting.ca","tsvader@meierinc.com","charles.menun@menunrisk.com","miguel@mgelcor.com","cwu@mhpse.com","mmusso@m-west-assoc.com","agilani@miyamotointernational.com","yilmazcem01@gmail.com","mkuleli@gmail.com","k_achyut@hotmail.com","kocakyt42@hotmail.com","khans2@mmm.ca","pfendtf@mmm.ca","caleb@morriseng.com","kvanluchene@m-m.net","prachakatla@melick-tully.com","peter.gunness@mts.com","justin.c.beutel@mwhglobal.com","christine.t.weber@mwhglobal.com","dina.b.hunt@mwhglobal.com","ali.amini@shaw.ca","ernestn@shaw.ca","rmiller@nautilusgrp.com","mtaylor2@dot.state.nv.us","jwieser@dot.state.nv.us","sed@newalbiongeotechnical.com","alan.witthoeft@gmail.com","dchu@ninyoandmoore.com","smarcinek@ninyoandmoore.com","raul.uribe@nist.gov","anne.hulsey@nist.gov","jeradhoffman@gmail.com","jhoffman@nwgeotech.com","brian.queen@epa.ohio.gov","gokturkonem@gmail.com","campbell.keepa@opus.co","adhikarig@gmail.com","aleksandarzurovski@yahoo.com","brsshn85@gmail.com","chenj2@pbworld.com","habibghodsi@yahoo.com","rymann@gmail.com","aschildmeyer@pcs-structural.com","pbcs911@gmail.com","matthew.dorsey@psiusa.com","rtsewell@aol.com","rajiv@rhinooneeng.com","cengizerdogan1061@hotmail.com","sidorvash@gmail.com","bmkim81@gmail.com","timothydancheta@gmail.com","emel.seyhan@gmail.com","byungmin.kim@rms.com","maxemow@silman.com","ldmedeiros@rockridgegeo.com","wmcvitty@ruthcheck.com","mschotanus@ruthchek.com","ljiang@ruthchek.com","lburkett@ruthchek.com","dmurphy@sageengineers.com","mryan@sageengineers.com","sukyong.moon@samsung.com","yk0127.park@samsung.com","asad.bassam@gmail.com","gunupkwon@yahoo.com","bkosbab@scsolutions.com","thastings@schnabel-eng.com","ozhang@schramminc.com","matthews@s-frame.com","setoodehz36@gmail.com","ashahbzn@gmail.com","laserlipz@gmail.com","mja@shanwil.com","kpc@shanwil.com","jzb@shanwil.com","wetzelnick@gmail.com","atsarawit@sgh.com","aanup@sgh.com","srshepherd@sgh.com","scdalring@sgh.com","rsoto@sgh.com","rohamburger@sgh.com","rcappa@sgh.com","ibaig@sgh.com","danielbernard@bell.net","bsoumountha@sbcglobal.net","matthew.muto@sce.com","ttnguyen@southernco.com","x2bsurfa@southernco.com","x2skanji@southernco.com","syang@srk.com","adam.mcintyre@stantec.com","chris.longley@stantec.com","april.welshans@stantec.com","cpapadelis@vecsa.com","tradford@vecsa.com","mwang@vecsa.com","johnf@sens-usa.com","gpmurphy@tectonicengineering.com","byoung@tectonicengineering.com","huffte@gmail.com","jsadler@terra-associates.com","dabaska@terracon.com","rwiddle@terrapower.com","duan.peng@tetratech.com","lalinda.weerasekara@tetratech.com","zach.wilder@tetratech.com","danieljmolinaro@gmail.com","ebusby7@gmail.com","jtobolski@gmail.com","drumfanatic@gmail.com","ckirchner@thorntontomasetti.com","cgulec@thorntontomasetti.com","dma@thurber.ca","cng@thurber.ca","mbosse@thurber.ca","iriveracruz@thurber.ca","maksim.yegorov@tippingmar.com","munsoncliff@gmail.com","jeff.bayless@urs.com","mark.dober@urs.com","mehrdad.hosseini@urs.com","andreas.skarlatoudis@urs.com","summitstructural@msn.com","pliu@usbr.gov","thomas.weaver@nrc.gov","scott.stovall@nrc.gov","josh.m.corbett@usace.army.mil","george.hu@usace.army.mil","stephenromero@fs.fed.us","gord461@ecy.wa.gov","tbajwa@westerngeo.ca","jmeng@westerngeo.ca","joel.aguilar@worleyparsons.com","sarahk@wrkengrs.com","palmerp@wsdot.wa.gov","saeedtowfigh@gmail.com","afyodorova@yu-associates.com","ethanm@zfa.com","caseyc@zfa.com"]

  #   verifiedEmailList.each do |userEmail|
  #     prepared_stmt = "UPDATE users SET users.userType = \"VerifiedUser\" WHERE (users.email=\"#{userEmail}\" and users.userType <> \"Administrator\");"
  #     logger.info prepared_stmt
  #       (ActiveRecord::Base.connection.execute(prepared_stmt))
  #       ActiveRecord::Base.connection_handler.clear_active_connections!
  #   end

  # end



  def AdminFunction_ManageBannedUser()
    bannedUserIDList = [172,173,174,175,176,179,180,181,182,183,184,185]

    bannedUserIDList.each do |bannedUserID|
      prepared_stmt = "UPDATE users SET users.download_limit_2week = 0 WHERE (users.id=#{bannedUserID});"
      (ActiveRecord::Base.connection.execute(prepared_stmt))
      prepared_stmt = "UPDATE users SET users.download_limit_1month = 0 WHERE (users.id=#{bannedUserID});"
      (ActiveRecord::Base.connection.execute(prepared_stmt))
      prepared_stmt = "UPDATE users SET users.userType = \"BannedUser\" WHERE (users.id=#{bannedUserID});"
      (ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!
    end

  end


  def AdminFunctions_WriteStats()


    puts 'in AdminFunctions_WriteStats'

    # startTime=Time.now
    unique_token=self.unique_token

    verifiedEmailList = self.imageName
    verifiedEmailList = verifiedEmailList.split(',')
    verifiedEmailList.each do |userEmail|
      prepared_stmt = "UPDATE users SET users.userType = \"VerifiedUser\" WHERE (users.email=\"#{userEmail}\" and users.userType <> \"Administrator\");"
      logger.info prepared_stmt
      (ActiveRecord::Base.connection.execute(prepared_stmt))
      ActiveRecord::Base.connection_handler.clear_active_connections!
    end


    outputLines = []

    #prepared_stmt ="delete z4_file_downloads.* from z4_file_downloads where ((z4_file_downloads.user_id =5 or z4_file_downloads.user_id =16)); "
    #   ActiveRecord::Base.connection.execute(prepared_stmt)
    #     ActiveRecord::Base.connection_handler.clear_active_connections!


    #prepared_stmt ="delete spectras.* from spectras where ((spectras.user_id =5 or spectras.user_id =16)); "
    #   ActiveRecord::Base.connection.execute(prepared_stmt)
    #     ActiveRecord::Base.connection_handler.clear_active_connections!

    #prepared_stmt ="delete searches.* from searches where ((searches.user_id =5 or searches.user_id =16)); "
    #   ActiveRecord::Base.connection.execute(prepared_stmt)
    #     ActiveRecord::Base.connection_handler.clear_active_connections!
    # commandList = [
    #  # "update users set encrypted_password = 55 where id=4; ",
    #  # "update users set password_salt = 55 where id=4; ",
    #  "update users set download_limit_2week = null where id=4; ",
    #  "update users set download_limit_1month = null where id=4; ",
    #  "update users set userType = \"User\" where id=4; ",
    # ]

    # commandList.each do |prepared_stmt|
    #     (ActiveRecord::Base.connection.execute(prepared_stmt))
    # end
    #       ActiveRecord::Base.connection_handler.clear_active_connections!




    inList = [
        "",
        "  -- PEER Ground Motion Database Administrator Report -- NGA-West2 -- " + Date.today.to_s,
        "",
    ]
    inList.each do |prepared_stmt|
      outputLines << prepared_stmt
    end

    puts 10


    thisTable = "users"
    skipList = ["encrypted_password", "password_salt"]
    inList = [
        "",
        "  -- #{thisTable} -- " + Date.today.to_s,
        "",
    ]
    inList.each do |prepared_stmt|
      outputLines << prepared_stmt
    end
    prepared_stmt = "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='rdgml_development' AND `TABLE_NAME`='#{thisTable}' ORDER BY ORDINAL_POSITION;"
    result_setColumnName=(ActiveRecord::Base.connection.execute(prepared_stmt))
    ActiveRecord::Base.connection_handler.clear_active_connections!
    columnLabelString = ""
    skipCounterList = []
    thisCounter = 1
    result_setColumnName.each{|thisColumnLabel|
      if skipList.include?(thisColumnLabel[0]) == true
        skipCounterList << thisCounter
      else
        columnLabelString = "#{columnLabelString},#{thisColumnLabel[0]}"
      end
      thisCounter +=1
    }
    outputLines << columnLabelString
    prepared_stmt = "SELECT * FROM #{thisTable};"
    result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
    ActiveRecord::Base.connection_handler.clear_active_connections!
    result_set.each{|rowValue|
      thisCounter = 1
      thisRowString = ""
      rowValue.each{|colValue|
        if skipCounterList.include?(thisCounter) == false
          if (colValue.is_a? Array)
            colValue = colValue.join(" ")
          end
          if (colValue.is_a? String)
            colValue.strip!
            colValue = colValue.gsub(",", " ")
            colValue = colValue.gsub("\t", " ")
          end
          thisRowString = "#{thisRowString},#{colValue}"
        end
        thisCounter +=1
      }
      outputLines << thisRowString
    }
    puts 20

    thisTable = "spectras"
    skipList = []
    inList = [
        "",
        "  -- #{thisTable} --" + Date.today.to_s,
        "",
    ]
    inList.each do |prepared_stmt|
      outputLines << prepared_stmt
    end
    prepared_stmt = "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='rdgml_development' AND `TABLE_NAME`='#{thisTable}' ORDER BY ORDINAL_POSITION;"
    result_setColumnName=(ActiveRecord::Base.connection.execute(prepared_stmt))
    ActiveRecord::Base.connection_handler.clear_active_connections!
    columnLabelString = "user_id"
    skipCounterList = []
    thisCounter = 1
    user_idCounter = 1
    result_setColumnName.each{|thisColumnLabel|
      if skipList.include?(thisColumnLabel[0]) == true
        skipCounterList << thisCounter
      else
        columnLabelString = "#{columnLabelString},#{thisColumnLabel[0]}"
      end
      if thisColumnLabel[0] == "user_id"
        user_idCounter = thisCounter
      end
      thisCounter +=1
    }
    outputLines << columnLabelString
    prepared_stmt = "SELECT * FROM #{thisTable};"
    result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
    ActiveRecord::Base.connection_handler.clear_active_connections!
    result_set.each{|rowValue|
      thisCounter = 1
      thisRowString = "#{rowValue[user_idCounter-1]}"
      rowValue.each{|colValue|
        if skipCounterList.include?(thisCounter) == false
          if (colValue.is_a? Array)
            colValue = colValue.join(" ")
          end
          if (colValue.is_a? String)
            colValue.strip!
            colValue = colValue.gsub(",", " ")
            colValue = colValue.gsub("\t", " ")
          end
          thisRowString = "#{thisRowString},#{colValue}"
        end
        thisCounter +=1
      }
      outputLines << thisRowString
    }
    puts 30

    thisTable = "searches"
    skipList = ["spectraInputArray","searchInputArray","scaleFactorArray","mseArray"]

    inList = [
        "",
        "  -- #{thisTable} -- " + Date.today.to_s,
        "",
    ]

    inList.each do |prepared_stmt|
      outputLines << prepared_stmt
    end
    prepared_stmt = "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='rdgml_development' AND `TABLE_NAME`='#{thisTable}' ORDER BY ORDINAL_POSITION;"
    result_setColumnName=(ActiveRecord::Base.connection.execute(prepared_stmt))
    ActiveRecord::Base.connection_handler.clear_active_connections!
    columnLabelString = ""
    skipCounterList = []
    thisCounter = 1
    result_setColumnName.each{|thisColumnLabel|
      if skipList.include?(thisColumnLabel[0]) == true
        skipCounterList << thisCounter
      else
        columnLabelString = "#{columnLabelString},#{thisColumnLabel[0]}"
      end
      thisCounter +=1
    }
    outputLines << columnLabelString
    prepared_stmt = "SELECT * FROM #{thisTable};"
    result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
    ActiveRecord::Base.connection_handler.clear_active_connections!
    result_set.each{|rowValue|
      thisCounter = 1
      thisRowString = ""
      rowValue.each{|colValue|
        if skipCounterList.include?(thisCounter) == false
          if (colValue.is_a? Array)
            colValue = colValue.join(" ")
          end
          colValue = colValue.to_s
          colValue.strip!
          colValue = colValue.gsub(",", " ")
          colValue = colValue.gsub("\t", " ")
          thisRowString = "#{thisRowString},#{colValue}"
        end
        thisCounter +=1
      }
      outputLines << thisRowString
    }


    puts 40

    thisTable = "z4_file_downloads"
    skipList = []
    inList = [
        "",
        "  -- #{thisTable} -- " + Date.today.to_s,
        "",
    ]
    inList.each do |prepared_stmt|
      outputLines << prepared_stmt
    end
    prepared_stmt = "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='rdgml_development' AND `TABLE_NAME`='#{thisTable}' ORDER BY ORDINAL_POSITION;"
    result_setColumnName=(ActiveRecord::Base.connection.execute(prepared_stmt))
    ActiveRecord::Base.connection_handler.clear_active_connections!
    columnLabelString = ""
    skipCounterList = []
    thisCounter = 1
    result_setColumnName.each{|thisColumnLabel|
      if skipList.include?(thisColumnLabel[0]) == true
        skipCounterList << thisCounter
      else
        columnLabelString = "#{columnLabelString},#{thisColumnLabel[0]}"
      end
      thisCounter +=1
    }
    outputLines << columnLabelString
    prepared_stmt = "SELECT * FROM #{thisTable};"
    result_set=(ActiveRecord::Base.connection.execute(prepared_stmt))
    ActiveRecord::Base.connection_handler.clear_active_connections!
    result_set.each{|rowValue|
      thisCounter = 1
      thisRowString = ""
      rowValue.each{|colValue|
        if skipCounterList.include?(thisCounter) == false
          if (colValue.is_a? Array)
            colValue = colValue.join(" ")
          end
          if (colValue.is_a? String)
            colValue.strip!
            colValue = colValue.gsub(",", " ")
            colValue = colValue.gsub("\t", " ")
          end
          thisRowString = "#{thisRowString},#{colValue}"
        end
        thisCounter +=1
      }
      outputLines << thisRowString
    }

    thisFilename = "NGAWest2AdminStats.txt"

    # directory = "/archive/ngawest2_files/Utility"
    # directory = "#{ENV['ARCHIVE_FILES_PATH']}Utility"
    directory_download  = "#{ENV['ARCHIVE_FILES_DOWNLOAD_PATH']}Utility"

    thisFullPath = File.join(directory_download, thisFilename)
    thisFullPath = File.expand_path(thisFullPath)

    fout = File.open(thisFullPath, 'w')
    fout.puts(outputLines)
    fout.close

    return 1

  end

end
