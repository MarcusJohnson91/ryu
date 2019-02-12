// Copyright 2018 Ulf Adams
//
// The contents of this file may be used under the terms of the Apache License,
// Version 2.0.
//
//    (See accompanying file LICENSE-Apache or copy at
//     http://www.apache.org/licenses/LICENSE-2.0)
//
// Alternatively, the contents of this file may be used under the terms of
// the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE-Boost or copy at
//     https://www.boost.org/LICENSE_1_0.txt)
//
// Unless required by applicable law or agreed to in writing, this software
// is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.
#ifndef RYU_D2S_FULL_TABLE_H
#define RYU_D2S_FULL_TABLE_H

#include "common.h"

// These tables are generated by PrintDoubleLookupTable.
static const uint64_t DOUBLE_POW5_INV_SPLIT[292][2] = {
    {                    1u, 288230376151711744u }, {  3689348814741910324u, 230584300921369395u },
    {  2951479051793528259u, 184467440737095516u }, { 17118578500402463900u, 147573952589676412u },
    { 12632330341676300947u, 236118324143482260u }, { 10105864273341040758u, 188894659314785808u },
    { 15463389048156653253u, 151115727451828646u }, { 17362724847566824558u, 241785163922925834u },
    { 17579528692795369969u, 193428131138340667u }, {  6684925324752475329u, 154742504910672534u },
    { 18074578149087781173u, 247588007857076054u }, { 18149011334012135262u, 198070406285660843u },
    {  3451162622983977240u, 158456325028528675u }, {  5521860196774363583u, 253530120045645880u },
    {  4417488157419490867u, 202824096036516704u }, {  7223339340677503017u, 162259276829213363u },
    {  7867994130342094503u, 259614842926741381u }, {  2605046489531765280u, 207691874341393105u },
    {  2084037191625412224u, 166153499473114484u }, { 10713157136084480204u, 265845599156983174u },
    { 12259874523609494487u, 212676479325586539u }, { 13497248433629505913u, 170141183460469231u },
    { 14216899864323388813u, 272225893536750770u }, { 11373519891458711051u, 217780714829400616u },
    {  5409467098425058518u, 174224571863520493u }, {  4965798542738183305u, 278759314981632789u },
    {  7661987648932456967u, 223007451985306231u }, {  2440241304404055250u, 178405961588244985u },
    {  3904386087046488400u, 285449538541191976u }, { 17880904128604832013u, 228359630832953580u },
    { 14304723302883865611u, 182687704666362864u }, { 15133127457049002812u, 146150163733090291u },
    { 16834306301794583852u, 233840261972944466u }, {  9778096226693756759u, 187072209578355573u },
    { 15201174610838826053u, 149657767662684458u }, {  2185786488890659746u, 239452428260295134u },
    {  5437978005854438120u, 191561942608236107u }, { 15418428848909281466u, 153249554086588885u },
    {  6222742084545298729u, 245199286538542217u }, { 16046240111861969953u, 196159429230833773u },
    {  1768945645263844993u, 156927543384667019u }, { 10209010661905972635u, 251084069415467230u },
    {  8167208529524778108u, 200867255532373784u }, { 10223115638361732810u, 160693804425899027u },
    {  1599589762411131202u, 257110087081438444u }, {  4969020624670815285u, 205688069665150755u },
    {  3975216499736652228u, 164550455732120604u }, { 13739044029062464211u, 263280729171392966u },
    {  7301886408508061046u, 210624583337114373u }, { 13220206756290269483u, 168499666669691498u },
    { 17462981995322520850u, 269599466671506397u }, {  6591687966774196033u, 215679573337205118u },
    { 12652048002903177473u, 172543658669764094u }, {  9175230360419352987u, 276069853871622551u },
    {  3650835473593572067u, 220855883097298041u }, { 17678063637842498946u, 176684706477838432u },
    { 13527506561580357021u, 282695530364541492u }, {  3443307619780464970u, 226156424291633194u },
    {  6443994910566282300u, 180925139433306555u }, {  5155195928453025840u, 144740111546645244u },
    { 15627011115008661990u, 231584178474632390u }, { 12501608892006929592u, 185267342779705912u },
    {  2622589484121723027u, 148213874223764730u }, {  4196143174594756843u, 237142198758023568u },
    { 10735612169159626121u, 189713759006418854u }, { 12277838550069611220u, 151771007205135083u },
    { 15955192865369467629u, 242833611528216133u }, {  1696107848069843133u, 194266889222572907u },
    { 12424932722681605476u, 155413511378058325u }, {  1433148282581017146u, 248661618204893321u },
    { 15903913885032455010u, 198929294563914656u }, {  9033782293284053685u, 159143435651131725u },
    { 14454051669254485895u, 254629497041810760u }, { 11563241335403588716u, 203703597633448608u },
    { 16629290697806691620u, 162962878106758886u }, {   781423413297334329u, 260740604970814219u },
    {  4314487545379777786u, 208592483976651375u }, {  3451590036303822229u, 166873987181321100u },
    {  5522544058086115566u, 266998379490113760u }, {  4418035246468892453u, 213598703592091008u },
    { 10913125826658934609u, 170878962873672806u }, { 10082303693170474728u, 273406340597876490u },
    {  8065842954536379782u, 218725072478301192u }, { 17520720807854834795u, 174980057982640953u },
    {  5897060404116273733u, 279968092772225526u }, {  1028299508551108663u, 223974474217780421u },
    { 15580034865808528224u, 179179579374224336u }, { 17549358155809824511u, 286687326998758938u },
    {  2971440080422128639u, 229349861599007151u }, { 17134547323305344204u, 183479889279205720u },
    { 13707637858644275364u, 146783911423364576u }, { 14553522944347019935u, 234854258277383322u },
    {  4264120725993795302u, 187883406621906658u }, { 10789994210278856888u, 150306725297525326u },
    {  9885293106962350374u, 240490760476040522u }, {   529536856086059653u, 192392608380832418u },
    {  7802327114352668369u, 153914086704665934u }, {  1415676938738538420u, 246262538727465495u },
    {  1132541550990830736u, 197010030981972396u }, { 15663428499760305882u, 157608024785577916u },
    { 17682787970132668764u, 252172839656924666u }, { 10456881561364224688u, 201738271725539733u },
    { 15744202878575200397u, 161390617380431786u }, { 17812026976236499989u, 258224987808690858u },
    {  3181575136763469022u, 206579990246952687u }, { 13613306553636506187u, 165263992197562149u },
    { 10713244041592678929u, 264422387516099439u }, { 12259944048016053467u, 211537910012879551u },
    {  6118606423670932450u, 169230328010303641u }, {  2411072648389671274u, 270768524816485826u },
    { 16686253377679378312u, 216614819853188660u }, { 13349002702143502650u, 173291855882550928u },
    { 17669055508687693916u, 277266969412081485u }, { 14135244406950155133u, 221813575529665188u },
    {   240149081334393137u, 177450860423732151u }, { 11452284974360759988u, 283921376677971441u },
    {  5472479164746697667u, 227137101342377153u }, { 11756680961281178780u, 181709681073901722u },
    {  2026647139541122378u, 145367744859121378u }, { 18000030682233437097u, 232588391774594204u },
    { 18089373360528660001u, 186070713419675363u }, {  3403452244197197031u, 148856570735740291u },
    { 16513570034941246220u, 238170513177184465u }, { 13210856027952996976u, 190536410541747572u },
    {  3189987192878576934u, 152429128433398058u }, {  1414630693863812771u, 243886605493436893u },
    {  8510402184574870864u, 195109284394749514u }, { 10497670562401807014u, 156087427515799611u },
    {  9417575270359070576u, 249739884025279378u }, { 14912757845771077107u, 199791907220223502u },
    {  4551508647133041040u, 159833525776178802u }, { 10971762650154775986u, 255733641241886083u },
    { 16156107749607641435u, 204586912993508866u }, {  9235537384944202825u, 163669530394807093u },
    { 11087511001168814197u, 261871248631691349u }, { 12559357615676961681u, 209496998905353079u },
    { 13736834907283479668u, 167597599124282463u }, { 18289587036911657145u, 268156158598851941u },
    { 10942320814787415393u, 214524926879081553u }, { 16132554281313752961u, 171619941503265242u },
    { 11054691591134363444u, 274591906405224388u }, { 16222450902391311402u, 219673525124179510u },
    { 12977960721913049122u, 175738820099343608u }, { 17075388340318968271u, 281182112158949773u },
    {  2592264228029443648u, 224945689727159819u }, {  5763160197165465241u, 179956551781727855u },
    {  9221056315464744386u, 287930482850764568u }, { 14755542681855616155u, 230344386280611654u },
    { 15493782960226403247u, 184275509024489323u }, {  1326979923955391628u, 147420407219591459u },
    {  9501865507812447252u, 235872651551346334u }, { 11290841220991868125u, 188698121241077067u },
    {  1653975347309673853u, 150958496992861654u }, { 10025058185179298811u, 241533595188578646u },
    {  4330697733401528726u, 193226876150862917u }, { 14532604630946953951u, 154581500920690333u },
    {  1116074521063664381u, 247330401473104534u }, {  4582208431592841828u, 197864321178483627u },
    { 14733813189500004432u, 158291456942786901u }, { 16195403473716186445u, 253266331108459042u },
    {  5577625149489128510u, 202613064886767234u }, {  8151448934333213131u, 162090451909413787u },
    { 16731667109675051333u, 259344723055062059u }, { 17074682502481951390u, 207475778444049647u },
    {  6281048372501740465u, 165980622755239718u }, {  6360328581260874421u, 265568996408383549u },
    {  8777611679750609860u, 212455197126706839u }, { 10711438158542398211u, 169964157701365471u },
    {  9759603424184016492u, 271942652322184754u }, { 11497031554089123517u, 217554121857747803u },
    { 16576322872755119460u, 174043297486198242u }, { 11764721337440549842u, 278469275977917188u },
    { 16790474699436260520u, 222775420782333750u }, { 13432379759549008416u, 178220336625867000u },
    {  3045063541568861850u, 285152538601387201u }, { 17193446092222730773u, 228122030881109760u },
    { 13754756873778184618u, 182497624704887808u }, { 18382503128506368341u, 145998099763910246u },
    {  3586563302416817083u, 233596959622256395u }, {  2869250641933453667u, 186877567697805116u },
    { 17052795772514404226u, 149502054158244092u }, { 12527077977055405469u, 239203286653190548u },
    { 17400360011128145022u, 191362629322552438u }, {  2852241564676785048u, 153090103458041951u },
    { 15631632947708587046u, 244944165532867121u }, {  8815957543424959314u, 195955332426293697u },
    { 18120812478965698421u, 156764265941034957u }, { 14235904707377476180u, 250822825505655932u },
    {  4010026136418160298u, 200658260404524746u }, { 17965416168102169531u, 160526608323619796u },
    {  2919224165770098987u, 256842573317791675u }, {  2335379332616079190u, 205474058654233340u },
    {  1868303466092863352u, 164379246923386672u }, {  6678634360490491686u, 263006795077418675u },
    {  5342907488392393349u, 210405436061934940u }, {  4274325990713914679u, 168324348849547952u },
    { 10528270399884173809u, 269318958159276723u }, { 15801313949391159694u, 215455166527421378u },
    {  1573004715287196786u, 172364133221937103u }, { 17274202803427156150u, 275782613155099364u },
    { 17508711057483635243u, 220626090524079491u }, { 10317620031244997871u, 176500872419263593u },
    { 12818843235250086271u, 282401395870821749u }, { 13944423402941979340u, 225921116696657399u },
    { 14844887537095493795u, 180736893357325919u }, { 15565258844418305359u, 144589514685860735u },
    {  6457670077359736959u, 231343223497377177u }, { 16234182506113520537u, 185074578797901741u },
    {  9297997190148906106u, 148059663038321393u }, { 11187446689496339446u, 236895460861314229u },
    { 12639306166338981880u, 189516368689051383u }, { 17490142562555006151u, 151613094951241106u },
    {  2158786396894637579u, 242580951921985771u }, { 16484424376483351356u, 194064761537588616u },
    {  9498190686444770762u, 155251809230070893u }, { 11507756283569722895u, 248402894768113429u },
    { 12895553841597688639u, 198722315814490743u }, { 17695140702761971558u, 158977852651592594u },
    { 17244178680193423523u, 254364564242548151u }, { 10105994129412828495u, 203491651394038521u },
    {  4395446488788352473u, 162793321115230817u }, { 10722063196803274280u, 260469313784369307u },
    {  1198952927958798777u, 208375451027495446u }, { 15716557601334680315u, 166700360821996356u },
    { 17767794532651667857u, 266720577315194170u }, { 14214235626121334286u, 213376461852155336u },
    {  7682039686155157106u, 170701169481724269u }, {  1223217053622520399u, 273121871170758831u },
    { 15735968901865657612u, 218497496936607064u }, { 16278123936234436413u, 174797997549285651u },
    {   219556594781725998u, 279676796078857043u }, {  7554342905309201445u, 223741436863085634u },
    {  9732823138989271479u, 178993149490468507u }, {   815121763415193074u, 286389039184749612u },
    { 11720143854957885429u, 229111231347799689u }, { 13065463898708218666u, 183288985078239751u },
    {  6763022304224664610u, 146631188062591801u }, {  3442138057275642729u, 234609900900146882u },
    { 13821756890046245153u, 187687920720117505u }, { 11057405512036996122u, 150150336576094004u },
    {  6623802375033462826u, 240240538521750407u }, { 16367088344252501231u, 192192430817400325u },
    { 13093670675402000985u, 153753944653920260u }, {  2503129006933649959u, 246006311446272417u },
    { 13070549649772650937u, 196805049157017933u }, { 17835137349301941396u, 157444039325614346u },
    {  2710778055689733971u, 251910462920982955u }, {  2168622444551787177u, 201528370336786364u },
    {  5424246770383340065u, 161222696269429091u }, {  1300097203129523457u, 257956314031086546u },
    { 15797473021471260058u, 206365051224869236u }, {  8948629602435097724u, 165092040979895389u },
    {  3249760919670425388u, 264147265567832623u }, {  9978506365220160957u, 211317812454266098u },
    { 15361502721659949412u, 169054249963412878u }, {  2442311466204457120u, 270486799941460606u },
    { 16711244431931206989u, 216389439953168484u }, { 17058344360286875914u, 173111551962534787u },
    { 12535955717491360170u, 276978483140055660u }, { 10028764573993088136u, 221582786512044528u },
    { 15401709288678291155u, 177266229209635622u }, {  9885339602917624555u, 283625966735416996u },
    {  4218922867592189321u, 226900773388333597u }, { 14443184738299482427u, 181520618710666877u },
    {  4175850161155765295u, 145216494968533502u }, { 10370709072591134795u, 232346391949653603u },
    { 15675264887556728482u, 185877113559722882u }, {  5161514280561562140u, 148701690847778306u },
    {   879725219414678777u, 237922705356445290u }, {   703780175531743021u, 190338164285156232u },
    { 11631070584651125387u, 152270531428124985u }, {   162968861732249003u, 243632850284999977u },
    { 11198421533611530172u, 194906280227999981u }, {  5269388412147313814u, 155925024182399985u },
    {  8431021459435702103u, 249480038691839976u }, {  3055468352806651359u, 199584030953471981u },
    { 17201769941212962380u, 159667224762777584u }, { 16454785461715008838u, 255467559620444135u },
    { 13163828369372007071u, 204374047696355308u }, { 17909760324981426303u, 163499238157084246u },
    {  2830174816776909822u, 261598781051334795u }, {  2264139853421527858u, 209279024841067836u },
    { 16568707141704863579u, 167423219872854268u }, {  4373838538276319787u, 267877151796566830u },
    {  3499070830621055830u, 214301721437253464u }, {  6488605479238754987u, 171441377149802771u },
    {  3003071137298187333u, 274306203439684434u }, {  6091805724580460189u, 219444962751747547u },
    { 15941491023890099121u, 175555970201398037u }, { 10748990379256517301u, 280889552322236860u },
    {  8599192303405213841u, 224711641857789488u }, { 14258051472207991719u, 179769313486231590u }
};

static const uint64_t DOUBLE_POW5_SPLIT[326][2] = {
    {                    0u,  72057594037927936u }, {                    0u,  90071992547409920u },
    {                    0u, 112589990684262400u }, {                    0u, 140737488355328000u },
    {                    0u,  87960930222080000u }, {                    0u, 109951162777600000u },
    {                    0u, 137438953472000000u }, {                    0u,  85899345920000000u },
    {                    0u, 107374182400000000u }, {                    0u, 134217728000000000u },
    {                    0u,  83886080000000000u }, {                    0u, 104857600000000000u },
    {                    0u, 131072000000000000u }, {                    0u,  81920000000000000u },
    {                    0u, 102400000000000000u }, {                    0u, 128000000000000000u },
    {                    0u,  80000000000000000u }, {                    0u, 100000000000000000u },
    {                    0u, 125000000000000000u }, {                    0u,  78125000000000000u },
    {                    0u,  97656250000000000u }, {                    0u, 122070312500000000u },
    {                    0u,  76293945312500000u }, {                    0u,  95367431640625000u },
    {                    0u, 119209289550781250u }, {  4611686018427387904u,  74505805969238281u },
    { 10376293541461622784u,  93132257461547851u }, {  8358680908399640576u, 116415321826934814u },
    {   612489549322387456u,  72759576141834259u }, { 14600669991935148032u,  90949470177292823u },
    { 13639151471491547136u, 113686837721616029u }, {  3213881284082270208u, 142108547152020037u },
    {  4314518811765112832u,  88817841970012523u }, {   781462496279003136u, 111022302462515654u },
    { 10200200157203529728u, 138777878078144567u }, { 13292654125893287936u,  86736173798840354u },
    {  7392445620511834112u, 108420217248550443u }, {  4628871007212404736u, 135525271560688054u },
    { 16728102434789916672u,  84703294725430033u }, {  7075069988205232128u, 105879118406787542u },
    { 18067209522111315968u, 132348898008484427u }, {  8986162942105878528u,  82718061255302767u },
    {  6621017659204960256u, 103397576569128459u }, {  3664586055578812416u, 129246970711410574u },
    { 16125424340018921472u,  80779356694631608u }, {  1710036351314100224u, 100974195868289511u },
    { 15972603494424788992u, 126217744835361888u }, {  9982877184015493120u,  78886090522101180u },
    { 12478596480019366400u,  98607613152626475u }, { 10986559581596820096u, 123259516440783094u },
    {  2254913720070624656u,  77037197775489434u }, { 12042014186943056628u,  96296497219361792u },
    { 15052517733678820785u, 120370621524202240u }, {  9407823583549262990u,  75231638452626400u },
    { 11759779479436578738u,  94039548065783000u }, { 14699724349295723422u, 117549435082228750u },
    {  4575641699882439235u,  73468396926392969u }, { 10331238143280436948u,  91835496157991211u },
    {  8302361660673158281u, 114794370197489014u }, {  1154580038986672043u, 143492962746861268u },
    {  9944984561221445835u,  89683101716788292u }, { 12431230701526807293u, 112103877145985365u },
    {  1703980321626345405u, 140129846432481707u }, { 17205888765512323542u,  87581154020301066u },
    { 12283988920035628619u, 109476442525376333u }, {  1519928094762372062u, 136845553156720417u },
    { 12479170105294952299u,  85528470722950260u }, { 15598962631618690374u, 106910588403687825u },
    {  5663645234241199255u, 133638235504609782u }, { 17374836326682913246u,  83523897190381113u },
    {  7883487353071477846u, 104404871487976392u }, {  9854359191339347308u, 130506089359970490u },
    { 10770660513014479971u,  81566305849981556u }, { 13463325641268099964u, 101957882312476945u },
    {  2994098996302961243u, 127447352890596182u }, { 15706369927971514489u,  79654595556622613u },
    {  5797904354682229399u,  99568244445778267u }, {  2635694424925398845u, 124460305557222834u },
    {  6258995034005762182u,  77787690973264271u }, {  3212057774079814824u,  97234613716580339u },
    { 17850130272881932242u, 121543267145725423u }, { 18073860448192289507u,  75964541966078389u },
    {  8757267504958198172u,  94955677457597987u }, {  6334898362770359811u, 118694596821997484u },
    { 13182683513586250689u,  74184123013748427u }, { 11866668373555425458u,  92730153767185534u },
    {  5609963430089506015u, 115912692208981918u }, { 17341285199088104971u,  72445432630613698u },
    { 12453234462005355406u,  90556790788267123u }, { 10954857059079306353u, 113195988485333904u },
    { 13693571323849132942u, 141494985606667380u }, { 17781854114260483896u,  88434366004167112u },
    {  3780573569116053255u, 110542957505208891u }, {   114030942967678664u, 138178696881511114u },
    {  4682955357782187069u,  86361685550944446u }, { 15077066234082509644u, 107952106938680557u },
    {  5011274737320973344u, 134940133673350697u }, { 14661261756894078100u,  84337583545844185u },
    {  4491519140835433913u, 105421979432305232u }, {  5614398926044292391u, 131777474290381540u },
    { 12732371365632458552u,  82360921431488462u }, {  6692092170185797382u, 102951151789360578u },
    { 17588487249587022536u, 128688939736700722u }, { 15604490549419276989u,  80430587335437951u },
    { 14893927168346708332u, 100538234169297439u }, { 14005722942005997511u, 125672792711621799u },
    { 15671105866394830300u,  78545495444763624u }, {  1142138259283986260u,  98181869305954531u },
    { 15262730879387146537u, 122727336632443163u }, {  7233363790403272633u,  76704585395276977u },
    { 13653390756431478696u,  95880731744096221u }, {  3231680390257184658u, 119850914680120277u },
    {  4325643253124434363u,  74906821675075173u }, { 10018740084832930858u,  93633527093843966u },
    {  3300053069186387764u, 117041908867304958u }, { 15897591223523656064u,  73151193042065598u },
    { 10648616992549794273u,  91438991302581998u }, {  4087399203832467033u, 114298739128227498u },
    { 14332621041645359599u, 142873423910284372u }, { 18181260187883125557u,  89295889943927732u },
    {  4279831161144355331u, 111619862429909666u }, { 14573160988285219972u, 139524828037387082u },
    { 13719911636105650386u,  87203017523366926u }, {  7926517508277287175u, 109003771904208658u },
    {   684774848491833161u, 136254714880260823u }, {  7345513307948477581u,  85159196800163014u },
    { 18405263671790372785u, 106448996000203767u }, { 18394893571310578077u, 133061245000254709u },
    { 13802651491282805250u,  83163278125159193u }, {  3418256308821342851u, 103954097656448992u },
    {  4272820386026678563u, 129942622070561240u }, {  2670512741266674102u,  81214138794100775u },
    { 17173198981865506339u, 101517673492625968u }, {  3019754653622331308u, 126897091865782461u },
    {  4193189667727651020u,  79310682416114038u }, { 14464859121514339583u,  99138353020142547u },
    { 13469387883465536574u, 123922941275178184u }, {  8418367427165960359u,  77451838296986365u },
    { 15134645302384838353u,  96814797871232956u }, {   471562554271496325u, 121018497339041196u },
    {  9518098633274461011u,  75636560836900747u }, {  7285937273165688360u,  94545701046125934u },
    { 18330793628311886258u, 118182126307657417u }, {  4539216990053847055u,  73863828942285886u },
    { 14897393274422084627u,  92329786177857357u }, {  4786683537745442072u, 115412232722321697u },
    { 14520892257159371055u,  72132645451451060u }, { 18151115321449213818u,  90165806814313825u },
    {  8853836096529353561u, 112707258517892282u }, {  1843923083806916143u, 140884073147365353u },
    { 12681666973447792349u,  88052545717103345u }, {  2017025661527576725u, 110065682146379182u },
    { 11744654113764246714u, 137582102682973977u }, {   422879793461572340u,  85988814176858736u },
    {   528599741826965425u, 107486017721073420u }, {   660749677283706782u, 134357522151341775u },
    {  7330497575943398595u,  83973451344588609u }, { 13774807988356636147u, 104966814180735761u },
    {  3383451930163631472u, 131208517725919702u }, { 15949715511634433382u,  82005323578699813u },
    {  6102086334260878016u, 102506654473374767u }, {  3015921899398709616u, 128133318091718459u },
    { 18025852251620051174u,  80083323807324036u }, {  4085571240815512351u, 100104154759155046u },
    { 14330336087874166247u, 125130193448943807u }, { 15873989082562435760u,  78206370905589879u },
    { 15230800334775656796u,  97757963631987349u }, {  5203442363187407284u, 122197454539984187u },
    {   946308467778435600u,  76373409087490117u }, {  5794571603150432404u,  95466761359362646u },
    { 16466586540792816313u, 119333451699203307u }, {  7985773578781816244u,  74583407312002067u },
    {  5370530955049882401u,  93229259140002584u }, {  6713163693812353001u, 116536573925003230u },
    { 18030785363914884337u,  72835358703127018u }, { 13315109668038829614u,  91044198378908773u },
    {  2808829029766373305u, 113805247973635967u }, { 17346094342490130344u, 142256559967044958u },
    {  6229622945628943561u,  88910349979403099u }, {  3175342663608791547u, 111137937474253874u },
    { 13192550366365765242u, 138922421842817342u }, {  3633657960551215372u,  86826513651760839u },
    { 18377130505971182927u, 108533142064701048u }, {  4524669058754427043u, 135666427580876311u },
    {  9745447189362598758u,  84791517238047694u }, {  2958436949848472639u, 105989396547559618u },
    { 12921418224165366607u, 132486745684449522u }, { 12687572408530742033u,  82804216052780951u },
    { 11247779492236039638u, 103505270065976189u }, {   224666310012885835u, 129381587582470237u },
    {  2446259452971747599u,  80863492239043898u }, { 12281196353069460307u, 101079365298804872u },
    { 15351495441336825384u, 126349206623506090u }, { 14206370669262903769u,  78968254139691306u },
    {  8534591299723853903u,  98710317674614133u }, { 15279925143082205283u, 123387897093267666u },
    { 14161639232853766206u,  77117435683292291u }, { 13090363022639819853u,  96396794604115364u },
    { 16362953778299774816u, 120495993255144205u }, { 12532689120651053212u,  75309995784465128u },
    { 15665861400813816515u,  94137494730581410u }, { 10358954714162494836u, 117671868413226763u },
    {  4168503687137865320u,  73544917758266727u }, {   598943590494943747u,  91931147197833409u },
    {  5360365506546067587u, 114913933997291761u }, { 11312142901609972388u, 143642417496614701u },
    {  9375932322719926695u,  89776510935384188u }, { 11719915403399908368u, 112220638669230235u },
    { 10038208235822497557u, 140275798336537794u }, { 10885566165816448877u,  87672373960336121u },
    { 18218643725697949000u, 109590467450420151u }, { 18161618638695048346u, 136988084313025189u },
    { 13656854658398099168u,  85617552695640743u }, { 12459382304570236056u, 107021940869550929u },
    {  1739169825430631358u, 133777426086938662u }, { 14922039196176308311u,  83610891304336663u },
    { 14040862976792997485u, 104513614130420829u }, {  3716020665709083144u, 130642017663026037u },
    {  4628355925281870917u,  81651261039391273u }, { 10397130925029726550u, 102064076299239091u },
    {  8384727637859770284u, 127580095374048864u }, {  5240454773662356427u,  79737559608780540u },
    {  6550568467077945534u,  99671949510975675u }, {  3576524565420044014u, 124589936888719594u },
    {  6847013871814915412u,  77868710555449746u }, { 17782139376623420074u,  97335888194312182u },
    { 13004302183924499284u, 121669860242890228u }, { 17351060901807587860u,  76043662651806392u },
    {  3242082053549933210u,  95054578314757991u }, { 17887660622219580224u, 118818222893447488u },
    { 11179787888887237640u,  74261389308404680u }, { 13974734861109047050u,  92826736635505850u },
    {  8245046539531533005u, 116033420794382313u }, { 16682369133275677888u,  72520887996488945u },
    {  7017903361312433648u,  90651109995611182u }, { 17995751238495317868u, 113313887494513977u },
    {  8659630992836983623u, 141642359368142472u }, {  5412269370523114764u,  88526474605089045u },
    { 11377022731581281359u, 110658093256361306u }, {  4997906377621825891u, 138322616570451633u },
    { 14652906532082110942u,  86451635356532270u }, {  9092761128247862869u, 108064544195665338u },
    {  2142579373455052779u, 135080680244581673u }, { 12868327154477877747u,  84425425152863545u },
    {  2250350887815183471u, 105531781441079432u }, {  2812938609768979339u, 131914726801349290u },
    {  6369772649532999991u,  82446704250843306u }, { 17185587848771025797u, 103058380313554132u },
    {  3035240737254230630u, 128822975391942666u }, {  6508711479211282048u,  80514359619964166u },
    { 17359261385868878368u, 100642949524955207u }, { 17087390713908710056u, 125803686906194009u },
    {  3762090168551861929u,  78627304316371256u }, {  4702612710689827411u,  98284130395464070u },
    { 15101637925217060072u, 122855162994330087u }, { 16356052730901744401u,  76784476871456304u },
    {  1998321839917628885u,  95980596089320381u }, {  7109588318324424010u, 119975745111650476u },
    { 13666864735807540814u,  74984840694781547u }, { 12471894901332038114u,  93731050868476934u },
    {  6366496589810271835u, 117163813585596168u }, {  3979060368631419896u,  73227383490997605u },
    {  9585511479216662775u,  91534229363747006u }, {  2758517312166052660u, 114417786704683758u },
    { 12671518677062341634u, 143022233380854697u }, {  1002170145522881665u,  89388895863034186u },
    { 10476084718758377889u, 111736119828792732u }, { 13095105898447972362u, 139670149785990915u },
    {  5878598177316288774u,  87293843616244322u }, { 16571619758500136775u, 109117304520305402u },
    { 11491152661270395161u, 136396630650381753u }, {   264441385652915120u,  85247894156488596u },
    {   330551732066143900u, 106559867695610745u }, {  5024875683510067779u, 133199834619513431u },
    { 10058076329834874218u,  83249896637195894u }, {  3349223375438816964u, 104062370796494868u },
    {  4186529219298521205u, 130077963495618585u }, { 14145795808130045513u,  81298727184761615u },
    { 13070558741735168987u, 101623408980952019u }, { 11726512408741573330u, 127029261226190024u },
    {  7329070255463483331u,  79393288266368765u }, { 13773023837756742068u,  99241610332960956u },
    { 17216279797195927585u, 124052012916201195u }, {  8454331864033760789u,  77532508072625747u },
    {  5956228811614813082u,  96915635090782184u }, {  7445286014518516353u, 121144543863477730u },
    {  9264989777501460624u,  75715339914673581u }, { 16192923240304213684u,  94644174893341976u },
    {  1794409976670715490u, 118305218616677471u }, {  8039035263060279037u,  73940761635423419u },
    {  5437108060397960892u,  92425952044279274u }, { 16019757112352226923u, 115532440055349092u },
    {   788976158365366019u,  72207775034593183u }, { 14821278253238871236u,  90259718793241478u },
    {  9303225779693813237u, 112824648491551848u }, { 11629032224617266546u, 141030810614439810u },
    { 11879831158813179495u,  88144256634024881u }, {  1014730893234310657u, 110180320792531102u },
    { 10491785653397664129u, 137725400990663877u }, {  8863209042587234033u,  86078375619164923u },
    {  6467325284806654637u, 107597969523956154u }, { 17307528642863094104u, 134497461904945192u },
    { 10817205401789433815u,  84060913690590745u }, { 18133192770664180173u, 105076142113238431u },
    { 18054804944902837312u, 131345177641548039u }, { 18201782118205355176u,  82090736025967524u },
    {  4305483574047142354u, 102613420032459406u }, { 14605226504413703751u, 128266775040574257u },
    {  2210737537617482988u,  80166734400358911u }, { 16598479977304017447u, 100208418000448638u },
    { 11524727934775246001u, 125260522500560798u }, {  2591268940807140847u,  78287826562850499u },
    { 17074144231291089770u,  97859783203563123u }, { 16730994270686474309u, 122324729004453904u },
    { 10456871419179046443u,  76452955627783690u }, {  3847717237119032246u,  95566194534729613u },
    {  9421332564826178211u, 119457743168412016u }, {  5888332853016361382u,  74661089480257510u },
    { 16583788103125227536u,  93326361850321887u }, { 16118049110479146516u, 116657952312902359u },
    { 16991309721690548428u,  72911220195563974u }, { 12015765115258409727u,  91139025244454968u },
    { 15019706394073012159u, 113923781555568710u }, {  9551260955736489391u, 142404726944460888u },
    {  5969538097335305869u,  89002954340288055u }, {  2850236603241744433u, 111253692925360069u }
};

#endif // RYU_D2S_FULL_TABLE_H
