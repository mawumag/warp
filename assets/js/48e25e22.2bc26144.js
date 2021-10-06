"use strict";(self.webpackChunkwebsite_2=self.webpackChunkwebsite_2||[]).push([[301],{3905:function(e,t,i){i.d(t,{Zo:function(){return c},kt:function(){return h}});var a=i(7294);function n(e,t,i){return t in e?Object.defineProperty(e,t,{value:i,enumerable:!0,configurable:!0,writable:!0}):e[t]=i,e}function o(e,t){var i=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);t&&(a=a.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),i.push.apply(i,a)}return i}function r(e){for(var t=1;t<arguments.length;t++){var i=null!=arguments[t]?arguments[t]:{};t%2?o(Object(i),!0).forEach((function(t){n(e,t,i[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(i)):o(Object(i)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(i,t))}))}return e}function l(e,t){if(null==e)return{};var i,a,n=function(e,t){if(null==e)return{};var i,a,n={},o=Object.keys(e);for(a=0;a<o.length;a++)i=o[a],t.indexOf(i)>=0||(n[i]=e[i]);return n}(e,t);if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(e);for(a=0;a<o.length;a++)i=o[a],t.indexOf(i)>=0||Object.prototype.propertyIsEnumerable.call(e,i)&&(n[i]=e[i])}return n}var s=a.createContext({}),p=function(e){var t=a.useContext(s),i=t;return e&&(i="function"==typeof e?e(t):r(r({},t),e)),i},c=function(e){var t=p(e.components);return a.createElement(s.Provider,{value:t},e.children)},u={inlineCode:"code",wrapper:function(e){var t=e.children;return a.createElement(a.Fragment,{},t)}},d=a.forwardRef((function(e,t){var i=e.components,n=e.mdxType,o=e.originalType,s=e.parentName,c=l(e,["components","mdxType","originalType","parentName"]),d=p(i),h=n,m=d["".concat(s,".").concat(h)]||d[h]||u[h]||o;return i?a.createElement(m,r(r({ref:t},c),{},{components:i})):a.createElement(m,r({ref:t},c))}));function h(e,t){var i=arguments,n=t&&t.mdxType;if("string"==typeof e||n){var o=i.length,r=new Array(o);r[0]=d;var l={};for(var s in t)hasOwnProperty.call(t,s)&&(l[s]=t[s]);l.originalType=e,l.mdxType="string"==typeof e?e:n,r[1]=l;for(var p=2;p<o;p++)r[p]=i[p];return a.createElement.apply(null,r)}return a.createElement.apply(null,i)}d.displayName="MDXCreateElement"},7350:function(e,t,i){i.r(t),i.d(t,{frontMatter:function(){return r},contentTitle:function(){return l},metadata:function(){return s},toc:function(){return p},default:function(){return u}});var a=i(7462),n=i(3366),o=(i(7294),i(3905)),r={sidebar_position:1},l="Best Practices for Building Data Processing Pipelines",s={unversionedId:"About_WARP/BestPractices",id:"About_WARP/BestPractices",isDocsHomePage:!1,title:"Best Practices for Building Data Processing Pipelines",description:"WARP pipeline development is guided by the best practices detailed below. We describe each of these best practices to give insight as to why they are important and we provide examples to give you a sense of how to apply them.",source:"@site/docs/About_WARP/BestPractices.md",sourceDirName:"About_WARP",slug:"/About_WARP/BestPractices",permalink:"/warp/docs/About_WARP/BestPractices",editUrl:"https://github.com/broadinstitute/warp/edit/develop/website/docs/About_WARP/BestPractices.md",tags:[],version:"current",lastUpdatedBy:"kew24",lastUpdatedAt:1633482549,formattedLastUpdatedAt:"10/6/2021",sidebarPosition:1,frontMatter:{sidebar_position:1},sidebar:"docsSidebar",previous:{title:"Welcome to WARP",permalink:"/warp/docs/get-started"},next:{title:"Pipeline Requirements",permalink:"/warp/docs/About_WARP/PipelineRequirements"}},p=[{value:"Automation",id:"automation",children:[{value:"What is Automation?",id:"what-is-automation",children:[]},{value:"Why do we care about automation?",id:"why-do-we-care-about-automation",children:[]}]},{value:"Testability",id:"testability",children:[{value:"What is a testable pipeline?",id:"what-is-a-testable-pipeline",children:[]},{value:"Why do we care about testabilty?",id:"why-do-we-care-about-testabilty",children:[]}]},{value:"Portability",id:"portability",children:[{value:"What is pipeline portability?",id:"what-is-pipeline-portability",children:[]},{value:"Why do we care about portability?",id:"why-do-we-care-about-portability",children:[]},{value:"Configurability for running on different technical infrastructures.",id:"configurability-for-running-on-different-technical-infrastructures",children:[]},{value:"Separation between the environment, the execution of the pipeline, and the pipeline itself.",id:"separation-between-the-environment-the-execution-of-the-pipeline-and-the-pipeline-itself",children:[]}]},{value:"Scaling Characteristics",id:"scaling-characteristics",children:[{value:"What do we mean by scaling characteristics?",id:"what-do-we-mean-by-scaling-characteristics",children:[]},{value:"Why care about scalability?",id:"why-care-about-scalability",children:[]}]},{value:"Maintainability",id:"maintainability",children:[{value:"What is a maintainable pipeline?",id:"what-is-a-maintainable-pipeline",children:[]},{value:"Why do we care about maintainability?",id:"why-do-we-care-about-maintainability",children:[]},{value:"Readability",id:"readability",children:[]},{value:"Modularity",id:"modularity",children:[]},{value:"Leveraging Standards",id:"leveraging-standards",children:[]}]},{value:"Versioning",id:"versioning",children:[]},{value:"Licensing",id:"licensing",children:[{value:"What is licensing?",id:"what-is-licensing",children:[]},{value:"Why do we care about licensing?",id:"why-do-we-care-about-licensing",children:[]}]}],c={toc:p};function u(e){var t=e.components,i=(0,n.Z)(e,["components"]);return(0,o.kt)("wrapper",(0,a.Z)({},c,i,{components:t,mdxType:"MDXLayout"}),(0,o.kt)("h1",{id:"best-practices-for-building-data-processing-pipelines"},"Best Practices for Building Data Processing Pipelines"),(0,o.kt)("p",null,"WARP pipeline development is guided by the best practices detailed below. We describe each of these best practices to give insight as to why they are important and we provide examples to give you a sense of how to apply them."),(0,o.kt)("p",null,"Overall, the best pipelines should be:"),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},"automated"),(0,o.kt)("li",{parentName:"ul"},"easily testable"),(0,o.kt)("li",{parentName:"ul"},"portable"),(0,o.kt)("li",{parentName:"ul"},"scalable to their data"),(0,o.kt)("li",{parentName:"ul"},"easy to maintain")),(0,o.kt)("h2",{id:"automation"},"Automation"),(0,o.kt)("h3",{id:"what-is-automation"},"What is Automation?"),(0,o.kt)("p",null,"Automation refers to the ability of a pipeline to run, end-to-end, without human intervention."),(0,o.kt)("h3",{id:"why-do-we-care-about-automation"},"Why do we care about automation?"),(0,o.kt)("p",null,"Pipelines cannot scale to large amounts of data, or many runs, if manual steps must be performed within the pipeline. They also cannot be part of an automated system if they in fact are not automated. Manual steps will bottleneck your entire system and can require unmanageable operations. Moreover, manual steps performed by humans will vary, and will promote the production of data that can not be appropriately harmonized."),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Reduce parameterization to minimal inputs that do not vary for each input data."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Remove the need for parameters, replacing them with data-driven settings."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Offer defaults that are generally applicable for inputs that cannot be defined in a data-driven manner."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Offer the ability to check the status of pipeline runs."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Don\u2019t"),": Assume any file produced at any step of the pipeline is ok. Always check the status of underlying tools (Eg. check return codes)."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Don\u2019t"),": Keep output files produced by steps of the pipeline that errored; people will accidently use them if they exist. (Do keep logs for debugging.)"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Don\u2019t"),": Delete outputs from steps that passed when the full pipeline fails, keeping them enables you to pick up where you left off."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Don\u2019t"),": Use tools that are \u201cbuggy\u201d or fragile, find alternatives or improve the tools.")),(0,o.kt)("h2",{id:"testability"},"Testability"),(0,o.kt)("h3",{id:"what-is-a-testable-pipeline"},"What is a testable pipeline?"),(0,o.kt)("p",null,"A testable pipeline is one in which isolated sections or the full pipeline can checked for specified characteristics without modifying the pipeline\u2019s code. Testability requires the existence of appropriate data with which to run the test and a testing checklist that reflects a clear understanding of how the data will be used to evaluate the pipeline."),(0,o.kt)("h3",{id:"why-do-we-care-about-testabilty"},"Why do we care about testabilty?"),(0,o.kt)("p",null,"The availability of test data enables validation that the pipeline can produce the desired outcome. Formulation of a testing checklist allows the developer to clearly define the capabilities of the pipeline and the parameters of its use."),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Provide example test data with your pipeline/tool."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Provide the results of an execution of your pipeline/tool on the test data."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Refer to at least one real data set appropriate for your tool/pipeline with example output from an execution of your pipeline or tool."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Provide a checker tool."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Include an automated testing suite for pipeline tasks."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Include automated tests for the pipeline as an integrated unit (pipeline benchmarking)."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Don\u2019t"),": Use tools that do not have automated testing suites."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Don\u2019t"),": Write tests that assume a specific instance of data.")),(0,o.kt)("h2",{id:"portability"},"Portability"),(0,o.kt)("h3",{id:"what-is-pipeline-portability"},"What is pipeline portability?"),(0,o.kt)("p",null,"Pipeline portability refers to the ability of a pipeline to execute successfully on multiple technical architectures."),(0,o.kt)("h3",{id:"why-do-we-care-about-portability"},"Why do we care about portability?"),(0,o.kt)("p",null,(0,o.kt)("em",{parentName:"p"},"Science.")," Science is not science if results are not reproducible; the scientific method cannot occur without a repeatable experiment that can be modified. Data processing pipelines are an essential part of some scientific inquiry and where they are leveraged they should be repeatable to validate and extend scientific discovery."),(0,o.kt)("p",null,(0,o.kt)("em",{parentName:"p"},"Impact.")," Pipelines will have greatest impact when they can be leveraged in multiple environments. The more technical requirements for installing and running of a pipeline, the longer it will take for a researcher to have a usable running pipeline."),(0,o.kt)("p",null,(0,o.kt)("em",{parentName:"p"},"Maintainability.")," Over the long term, it is easier to maintain pipelines that can be run in multiple environments. Portability avoids being tied to specific infrastructure and enables ease of deployment to development environments."),(0,o.kt)("p",null,"To ensure that others will be able to use your pipeline, avoid building in assumptions about environments and infrastructures in which it will run."),(0,o.kt)("h3",{id:"configurability-for-running-on-different-technical-infrastructures"},"Configurability for running on different technical infrastructures."),(0,o.kt)("p",null,"Code should not change to enable a pipeline to run on a different technical architecture; this change in execution environment should be configurable outside of the pipeline code."),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Use a workflow language that allows a separation between the code that executes the data processing logic and the logic to run the pipeline on an infrastructure. ",(0,o.kt)("a",{parentName:"li",href:"https://software.broadinstitute.org/wdl/documentation"},"WDL")," and ",(0,o.kt)("a",{parentName:"li",href:"https://www.commonwl.org/user_guide/rec-practices"},"CWL")," languages have this feature."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Don\u2019t"),": Put logic to run the pipeline in the same code that executes the logic to process the data.")),(0,o.kt)("h3",{id:"separation-between-the-environment-the-execution-of-the-pipeline-and-the-pipeline-itself"},"Separation between the environment, the execution of the pipeline, and the pipeline itself."),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Use a containerization technology, such as ",(0,o.kt)("a",{parentName:"li",href:"https://www.docker.com/"},"Docker"),", to execute software."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Incorporate into your testing the execution of the pipeline in multiple execution environments."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Don\u2019t"),": Put environmental paths in software tools or the workflow language. When they must exist they belong in the configuration or (if they refer to the execution environment) in the container\u2019s build instructions (Eg. Dockerfile).")),(0,o.kt)("h2",{id:"scaling-characteristics"},"Scaling Characteristics"),(0,o.kt)("h3",{id:"what-do-we-mean-by-scaling-characteristics"},"What do we mean by scaling characteristics?"),(0,o.kt)("p",null,"Scaling characteristics describe the performance of the pipeline given a certain amount of data. This is often described with Big O notation when describing algorithms. This answers the question: As the size of the data for the pipeline increases, how many additional computes are needed to process that data? One would want to avoid algorithms or tools that scale poorly, or improve this relationship to be linear (or better)."),(0,o.kt)("h3",{id:"why-care-about-scalability"},"Why care about scalability?"),(0,o.kt)("p",null,"If you have poor scaling characteristics, it may take an exponential amount of time to process more data. This will eventually require unreasonable amounts of time (and money if running in the cloud) and generally reduce the applicability of the pipeline."),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Measure the relationship between the size of your input (or something equivalent) and resources needed to successfully complete the pipeline."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Don\u2019t"),": Incorporate tools that have poor scaling characteristics unless they offer significant improvement to the pipeline.")),(0,o.kt)("h2",{id:"maintainability"},"Maintainability"),(0,o.kt)("h3",{id:"what-is-a-maintainable-pipeline"},"What is a maintainable pipeline?"),(0,o.kt)("p",null,"A pipeline that can be easily operated and updated is maintainable."),(0,o.kt)("h3",{id:"why-do-we-care-about-maintainability"},"Why do we care about maintainability?"),(0,o.kt)("p",null,"The majority of the life of code involves maintenance and updates. Design and initial implementation require vastly shorter amounts of time compared to the typical time period over which the code is operated and updated. This is generally true in many areas of software engineering. Thus it is important to engineer software so that the maintenance phase is manageable and does not burden new software development or operations."),(0,o.kt)("h3",{id:"readability"},"Readability"),(0,o.kt)("p",null,"Software is a living document that should be easily read and understood, regardless of who is the reader or author of the code."),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Work in a space where the code is easy to access and navigate (Eg. ",(0,o.kt)("a",{parentName:"li",href:"https://github.com/"},"GitHub"),")"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Use common software package structure and idioms to aid the navigation of the software."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Use automated documentation for all technical documents as much as possible."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Don\u2019t"),": Write a large amount of documentation that does not live beside or within the code itself (it will become out of date).")),(0,o.kt)("h3",{id:"modularity"},"Modularity"),(0,o.kt)("p",null,"Modularity enables small units of code to be independently benchmarked, updated, validated, and exchanged as science or technology changes. Using these small units enables more rapid updates and better adaptation to innovation."),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Save progress by creating intermediate output between modules as they successfully complete."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Don\u2019t"),": Make monolithic tasks that perform many functionalities for the sake of speed."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Don\u2019t"),": Break every functionality of a pipeline into a separate module. (This contrasts with not making monolithic tasks; there is an optimum between monolithic tasks and highly resolved modularity that is the goal. One can use benchmarking, the tendency for functionality to be updated, and how dependent functionalities are to get a sense of what should be separate and what can be combined.)")),(0,o.kt)("h3",{id:"leveraging-standards"},"Leveraging Standards"),(0,o.kt)("p",null,"We recommend using standard file formats and interfaces. In computational biology, ",(0,o.kt)("a",{parentName:"p",href:"https://www.ga4gh.org/genomic-data-toolkit/"},"GA4GH")," is a great source of these standards. In cases where new formats are needed, we recommend working with a standards group like ",(0,o.kt)("a",{parentName:"p",href:"https://www.ga4gh.org/"},"GA4GH")," if possible."),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": When using containerization technologies, follow best practices to assure associated images do not update without explicit updates."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Make both the images and the build files (Dockerfile) available to document the environment. ",(0,o.kt)("a",{parentName:"li",href:"https://docs.docker.com/develop/develop-images/dockerfile_best-practices/"},"More on Dockerfiles"),".")),(0,o.kt)("h2",{id:"versioning"},"Versioning"),(0,o.kt)("p",null,"Versioning pipelines and associated Docker images allows you to determine when and how data is created (provenance). As you make improvements and changes to your pipeline, it is important to know which version of the pipeline and software you used to create a given dataset so that it can be easily reproduced. This not only facilitates scientific reproducibility for the greater community, it also allows you to verify that new pipeline changes produce consistent results. We recommend choosing a consistent versioning system (for example, the ",(0,o.kt)("a",{parentName:"p",href:"https://semver.org/"},"semantic system"),") and tracking pipeline changes in a ",(0,o.kt)("a",{parentName:"p",href:"https://keepachangelog.com/en/1.0.0/"},"changelog"),"."),(0,o.kt)("h2",{id:"licensing"},"Licensing"),(0,o.kt)("h3",{id:"what-is-licensing"},"What is licensing?"),(0,o.kt)("p",null,'According to Wikipedia "A software license is a legal instrument (usually by way of contract law, with or without printed material) governing the use or redistribution of software.\u201d (see ',(0,o.kt)("a",{parentName:"p",href:"https://en.wikipedia.org/wiki/Software_license"},"this Wikipedia article")," for details)."),(0,o.kt)("h3",{id:"why-do-we-care-about-licensing"},"Why do we care about licensing?"),(0,o.kt)("div",{className:"admonition admonition-warning alert alert--danger"},(0,o.kt)("div",{parentName:"div",className:"admonition-heading"},(0,o.kt)("h5",{parentName:"div"},(0,o.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,o.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"12",height:"16",viewBox:"0 0 12 16"},(0,o.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M5.05.31c.81 2.17.41 3.38-.52 4.31C3.55 5.67 1.98 6.45.9 7.98c-1.45 2.05-1.7 6.53 3.53 7.7-2.2-1.16-2.67-4.52-.3-6.61-.61 2.03.53 3.33 1.94 2.86 1.39-.47 2.3.53 2.27 1.67-.02.78-.31 1.44-1.13 1.81 3.42-.59 4.78-3.42 4.78-5.56 0-2.84-2.53-3.22-1.25-5.61-1.52.13-2.03 1.13-1.89 2.75.09 1.08-1.02 1.8-1.86 1.33-.67-.41-.66-1.19-.06-1.78C8.18 5.31 8.68 2.45 5.05.32L5.03.3l.02.01z"}))),"NOTE")),(0,o.kt)("div",{parentName:"div",className:"admonition-content"},(0,o.kt)("p",{parentName:"div"},"this section is opinion and is NOT legal advice."))),(0,o.kt)("p",null,"Licenses sometimes legally bind you as to how you use tools, and sometimes the terms of the license transfer to the software and data that is produced. This can restrict the potential for leveraging the pipeline and may require additional work."),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Do"),": Select tools that are openly licensed to run in your pipelines to avoid the possibility that legal requirements will restrict execution where technical requirements do not."),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("em",{parentName:"li"},"Don\u2019t"),": Create software tools or libraries without licenses, clear guidance on your intent for use is important.")))}u.isMDXComponent=!0}}]);