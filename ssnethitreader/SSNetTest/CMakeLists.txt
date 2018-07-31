
simple_plugin(SSNetTest "module"
			larcore_Geometry_Geometry_service
                        larcore_Geometry
                        nusimdata_SimulationBase
			#larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        lardataobj_RecoBase
                        lardataobj_AnalysisBase 
                        #lardata_RecoObjects larreco_Deprecated  
                        #larreco_RecoAlg
                        #lardata_Utilities
			#lareventdisplay_EventDisplay
			#${LARRECO_LIB}
                        #${LARDATA_LIB}
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			#${LAREVENTDISPLAY_LIB}
			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
              
)

install_headers()
install_fhicl()
install_source( EXTRAS ADDITIONAL_NOTES.txt README.txt )
install_scripts()
