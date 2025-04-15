use serde_json::{Map, Value};
use std::any::Any;

/// Module descriptors
///
/// This module contains the definitions for
/// various module descriptors used in isotools.
/// Each descriptor represents a specific module
/// and contains various fields that can be set
/// and retrieved.
///
/// # Examples
///
/// ```rust, no_run
/// use isotools::modules::{ModuleDescriptor, ModuleType};
///
/// let module = ModuleDescriptor::with_schema(ModuleType::IntronRetention);
/// assert_eq!(module.get_value(Box::new(IntronRetentionValue::IsIntronRetention)), Some(Value::Bool(false)));
///
/// let mut module = module;
/// module.set_value(Box::new(IntronRetentionValue::IsIntronRetention), Value::Bool(true)).unwrap();
/// assert_eq!(module.get_value(Box::new(IntronRetentionValue::IsIntronRetention)), Some(Value::Bool(true)));
/// ```
#[derive(Debug)]
pub enum ModuleType {
    IntronRetention,
    StartTruncation,
    FusionDetection,
    PolyAPrediction,
}

/// ModuleMap trait
///
/// This trait defines the interface for module descriptors.
/// It allows for getting and setting values associated with
/// the module, as well as downcasting to specific types.
///
/// # Examples
///
/// ```rust, no_run
/// use isotools::modules::{ModuleMap, IntronRetentionDescriptor, IntronRetentionValue};
///
/// let mut module = IntronRetentionDescriptor::new();
/// module.set_value(Box::new(IntronRetentionValue::IsIntronRetention), Value::Bool(true)).unwrap();
/// assert_eq!(module.get_value(Box::new(IntronRetentionValue::IsIntronRetention)), Some(Value::Bool(true)));
///
/// let value = module.get_value(Box::new(IntronRetentionValue::IsIntronRetention));
/// assert_eq!(value, Some(Value::Bool(true)));
///
/// let value = module.get_value(Box::new(IntronRetentionValue::IsRetentionSupported));
/// assert_eq!(value, Some(Value::Null));
/// ```
pub trait ModuleMap: Any + Send + Sync {
    fn get_value(&self, key: Box<dyn Any>) -> Option<serde_json::Value>;
    fn set_value(&mut self, key: Box<dyn Any>, value: serde_json::Value) -> Result<(), String>;
    fn as_any(&self) -> &dyn Any;
    fn to_json(&self) -> Value;
}

/// Downcasting macro
///
/// This macro is used to downcast the module descriptor
/// to a specific type and format it for debugging.
///
/// # Examples
///
/// ```rust, no_run
/// use isotools::modules::{ModuleMap, IntronRetentionDescriptor};
///
/// let module = IntronRetentionDescriptor::new();
/// let result = downcast_dbg!(module, IntronRetentionDescriptor);
/// assert_eq!(result, "IntronRetentionDescriptor { ... }");
/// ```
macro_rules! downcast_dbg {
    ($formatter:expr, $module:expr, $($type:ty),+) => {
        {
            let mut result: Option<std::fmt::Result> = None;
            $(
                if result.is_none() {
                    if let Some(debuggable) = $module.as_any().downcast_ref::<$type>() {
                        result = Some(write!($formatter, "{:?}", debuggable));
                    }
                }
            )+
            result.unwrap_or_else(|| write!($formatter, "Unknown ModuleMap implementation"))
        }
    };
}

/// ModuleMap trait implementation for debugging
///
/// This implementation allows for formatting the module
/// descriptor for debugging purposes.
///
/// # Examples
///
/// ```rust, no_run
/// use isotools::modules::{ModuleMap, IntronRetentionDescriptor};
///
/// let module = IntronRetentionDescriptor::new();
/// let result = format!("{:?}", module);
/// assert_eq!(result, "IntronRetentionDescriptor { ... }");
/// ```
impl std::fmt::Debug for dyn ModuleMap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        downcast_dbg!(
            f,
            self,
            IntronRetentionDescriptor,
            StartTruncationDescriptor,
            FusionDetectionDescriptor,
            PolyAPredictionDescriptor
        )
    }
}

/// ModuleDescriptor struct
///
/// This struct represents a module descriptor
/// and contains a field for the module type.
///
/// # Examples
///
/// ```rust, no_run
/// use isotools::modules::{ModuleDescriptor, ModuleType};
///
/// let module = ModuleDescriptor::with_schema(ModuleType::IntronRetention);
/// assert_eq!(module.get_value(Box::new(IntronRetentionValue::IsIntronRetention)), Some(Value::Bool(false)));
/// ```
#[allow(dead_code)]
#[derive(Debug)]
pub struct ModuleDescriptor {
    module: ModuleType,
}

/// ModuleDescriptor implementation
///
/// This implementation provides a method to create
/// a new module descriptor with a specific schema.
impl ModuleDescriptor {
    /// Creates a new module descriptor with the specified schema.
    ///
    /// # Arguments
    ///
    /// * `module` - The module type to create a descriptor for.
    ///
    /// # Returns
    ///
    /// * A boxed module descriptor implementing the `ModuleMap` trait.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::{ModuleDescriptor, ModuleType};
    ///
    /// let module = ModuleDescriptor::with_schema(ModuleType::IntronRetention);
    /// assert_eq!(module.get_value(Box::new(IntronRetentionValue::IsIntronRetention)), Some(Value::Bool(false)));
    /// ```
    pub fn with_schema(module: ModuleType) -> Box<dyn ModuleMap> {
        match module {
            ModuleType::IntronRetention => IntronRetentionDescriptor::new(),
            ModuleType::StartTruncation => StartTruncationDescriptor::new(),
            ModuleType::FusionDetection => FusionDetectionDescriptor::new(),
            ModuleType::PolyAPrediction => PolyAPredictionDescriptor::new(),
        }
    }
}

/// IntronRetentionDescriptor struct
///
/// This struct represents the descriptor for
/// intron retention detection.
pub struct IntronRetentionDescriptor {
    pub intron_retention: Value,
    pub retention_support_type: Value,
    pub number_of_retentions: Value,
    pub coords_of_retention: Value,
    pub location_of_retention: Value,
    pub is_intron_retained_in_frame: Value,
    pub retains_rt_intron: Value,
    pub has_rt_intron: Value,
    pub has_rt_intron_map: Value,
    pub retention_acceptor_score: Value,
    pub retention_donor_score: Value,
    pub ref_introns_component_size: Value,
    pub query_component_size: Value,
    pub component_retention_ratio: Value,
    pub is_dirty_component: Value,
    pub exonic_status: Value,
    pub intronic_status: Value,
}

impl IntronRetentionDescriptor {
    /// Creates a new instance of the IntronRetentionDescriptor.
    ///
    /// # Returns
    ///
    /// * A boxed instance of the IntronRetentionDescriptor.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::IntronRetentionDescriptor;
    ///
    /// let descriptor = IntronRetentionDescriptor::new();
    /// assert_eq!(descriptor.intron_retention, Value::Bool(false));
    /// ```
    pub fn new() -> Box<Self> {
        Box::new(Self {
            intron_retention: Value::Bool(false),
            retention_support_type: Value::Null,
            number_of_retentions: Value::Number(0.into()),
            coords_of_retention: Value::Null,
            location_of_retention: Value::Null,
            is_intron_retained_in_frame: Value::Null,
            retains_rt_intron: Value::Null,
            has_rt_intron: Value::Null,
            has_rt_intron_map: Value::Null,
            retention_acceptor_score: Value::Null,
            retention_donor_score: Value::Null,
            ref_introns_component_size: Value::Null,
            query_component_size: Value::Null,
            component_retention_ratio: Value::Null,
            is_dirty_component: Value::Bool(false),
            exonic_status: Value::Null,
            intronic_status: Value::Null,
        })
    }
}

/// IntronRetentionValue enum
///
/// This enum defines the keys used to access
/// the values in the IntronRetentionDescriptor.
#[derive(Debug, Clone)]
pub enum IntronRetentionValue {
    IsIntronRetention,
    RetentionSupportType,
    NumberOfRetentions,
    RetentionCoords,
    RetentionLocation,
    IsIntronRetainedInFrame,
    RetainsRtIntron,
    HasRTIntron,
    HasRTIntronMap,
    RetentionAcceptorScore,
    RetentionDonorScore,
    RefIntronsComponentSize,
    QueryComponentSize,
    ComponentRetentionRatio,
    IsDirtyComponent,
    ExonicStatus,
    IntronicStatus,
}

/// ModuleMap trait implementation for IntronRetentionDescriptor
///
/// This implementation provides methods to get and set
/// values associated with the intron retention descriptor.
impl ModuleMap for IntronRetentionDescriptor {
    /// Gets the value associated with the specified key.
    ///
    /// # Arguments
    ///
    /// * `key` - The key to retrieve the value for.
    ///
    /// # Returns
    ///
    /// * An `Option` containing the value associated with the key,
    ///   or `None` if the key is not found.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::{ModuleMap, IntronRetentionDescriptor, IntronRetentionValue};
    ///
    /// let module = IntronRetentionDescriptor::new();
    /// assert_eq!(module.get_value(Box::new(IntronRetentionValue::IsIntronRetention)), Some(Value::Bool(false)));
    /// ```
    fn get_value(&self, key: Box<dyn Any>) -> Option<serde_json::Value> {
        if let Ok(key) = key.downcast::<IntronRetentionValue>() {
            match *key {
                IntronRetentionValue::IsIntronRetention => Some(self.intron_retention.clone()),
                IntronRetentionValue::RetentionSupportType => {
                    Some(self.retention_support_type.clone())
                }
                IntronRetentionValue::NumberOfRetentions => Some(self.number_of_retentions.clone()),
                IntronRetentionValue::RetentionCoords => Some(self.coords_of_retention.clone()),
                IntronRetentionValue::RetentionLocation => Some(self.location_of_retention.clone()),
                IntronRetentionValue::IsIntronRetainedInFrame => {
                    Some(self.is_intron_retained_in_frame.clone())
                }
                IntronRetentionValue::RetainsRtIntron => Some(self.retains_rt_intron.clone()),
                IntronRetentionValue::HasRTIntron => Some(self.has_rt_intron.clone()),
                IntronRetentionValue::HasRTIntronMap => Some(self.has_rt_intron_map.clone()),
                IntronRetentionValue::RetentionAcceptorScore => {
                    Some(self.retention_acceptor_score.clone())
                }
                IntronRetentionValue::RetentionDonorScore => {
                    Some(self.retention_donor_score.clone())
                }
                IntronRetentionValue::RefIntronsComponentSize => {
                    Some(self.ref_introns_component_size.clone())
                }
                IntronRetentionValue::QueryComponentSize => Some(self.query_component_size.clone()),
                IntronRetentionValue::ComponentRetentionRatio => {
                    Some(self.component_retention_ratio.clone())
                }
                IntronRetentionValue::IsDirtyComponent => Some(self.is_dirty_component.clone()),
                IntronRetentionValue::ExonicStatus => Some(self.exonic_status.clone()),
                IntronRetentionValue::IntronicStatus => Some(self.intronic_status.clone()),
            }
        } else {
            None
        }
    }

    /// Sets the value associated with the specified key.
    ///
    /// # Arguments
    ///
    /// * `key` - The key to set the value for.
    /// * `value` - The value to set.
    ///
    /// # Returns
    ///
    /// * A `Result` indicating success or failure.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::{ModuleMap, IntronRetentionDescriptor, IntronRetentionValue};
    ///
    /// let mut module = IntronRetentionDescriptor::new();
    /// module.set_value(Box::new(IntronRetentionValue::IsIntronRetention), Value::Bool(true)).unwrap();
    /// assert_eq!(module.get_value(Box::new(IntronRetentionValue::IsIntronRetention)), Some(Value::Bool(true)));
    /// ```
    #[inline(always)]
    fn set_value(&mut self, key: Box<dyn Any>, value: Value) -> Result<(), String> {
        if let Ok(key) = key.downcast::<IntronRetentionValue>() {
            match *key {
                IntronRetentionValue::IsIntronRetention => {
                    self.intron_retention = value;
                    Ok(())
                }
                IntronRetentionValue::RetentionSupportType => {
                    self.retention_support_type = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfRetentions => {
                    self.number_of_retentions = value;
                    Ok(())
                }
                IntronRetentionValue::RetentionCoords => {
                    self.coords_of_retention = value;
                    Ok(())
                }
                IntronRetentionValue::RetentionLocation => {
                    self.location_of_retention = value;
                    Ok(())
                }
                IntronRetentionValue::IsIntronRetainedInFrame => {
                    self.is_intron_retained_in_frame = value;
                    Ok(())
                }
                IntronRetentionValue::RetainsRtIntron => {
                    self.retains_rt_intron = value;
                    Ok(())
                }
                IntronRetentionValue::HasRTIntron => {
                    self.has_rt_intron = value;
                    Ok(())
                }
                IntronRetentionValue::HasRTIntronMap => {
                    self.has_rt_intron_map = value;
                    Ok(())
                }
                IntronRetentionValue::RetentionAcceptorScore => {
                    self.retention_acceptor_score = value;
                    Ok(())
                }
                IntronRetentionValue::RetentionDonorScore => {
                    self.retention_donor_score = value;
                    Ok(())
                }
                IntronRetentionValue::RefIntronsComponentSize => {
                    self.ref_introns_component_size = value;
                    Ok(())
                }
                IntronRetentionValue::QueryComponentSize => {
                    self.query_component_size = value;
                    Ok(())
                }
                IntronRetentionValue::ComponentRetentionRatio => {
                    self.component_retention_ratio = value;
                    Ok(())
                }
                IntronRetentionValue::IsDirtyComponent => {
                    self.is_dirty_component = value;
                    Ok(())
                }
                IntronRetentionValue::ExonicStatus => {
                    self.exonic_status = value;
                    Ok(())
                }
                IntronRetentionValue::IntronicStatus => {
                    self.intronic_status = value;
                    Ok(())
                }
            }
        } else {
            let err = format!("ERROR: You have tried to set a value for an unknown key!");
            log::error!("{}", err);
            Err(err)
        }
    }

    /// Downcasts the module descriptor to a specific type.
    ///
    /// # Returns
    ///
    /// * A reference to the module descriptor as a trait object.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::{ModuleMap, IntronRetentionDescriptor};
    ///
    /// let module = IntronRetentionDescriptor::new();
    /// assert_eq!(module.as_any().is::<IntronRetentionDescriptor>(), true);
    /// ```
    fn as_any(&self) -> &dyn Any {
        self
    }

    /// Converts the module descriptor to JSON format.
    ///
    /// # Returns
    ///
    /// * A JSON object representing the module descriptor.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::{ModuleMap, IntronRetentionDescriptor};
    ///
    /// let module = IntronRetentionDescriptor::new();
    /// let json = module.to_json();
    /// assert_eq!(json["is_intron_retention"], Value::Bool(false));
    /// ```
    fn to_json(&self) -> Value {
        let mut map = Map::new();

        macro_rules! insert {
            ($field:ident) => {
                map.insert(stringify!($field).to_string(), self.$field.clone());
            };
        }

        insert!(intron_retention);
        insert!(retention_support_type);
        insert!(number_of_retentions);
        insert!(coords_of_retention);
        insert!(location_of_retention);
        insert!(retains_rt_intron);
        insert!(has_rt_intron);
        insert!(has_rt_intron_map);
        insert!(is_intron_retained_in_frame);
        insert!(retention_acceptor_score);
        insert!(retention_donor_score);
        insert!(ref_introns_component_size);
        insert!(query_component_size);
        insert!(component_retention_ratio);
        insert!(is_dirty_component);
        insert!(exonic_status);
        insert!(intronic_status);

        Value::Object(map)
    }
}

/// Debug implementation for IntronRetentionDescriptor
///
/// This implementation provides a way to format the
/// IntronRetentionDescriptor for debugging purposes.
impl std::fmt::Debug for IntronRetentionDescriptor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{{
            is_intron_retention: {:?},
            retention_support_type: {:?},
            number_of_retentions: {:?},
            coords_of_retention: {:?},
            location_of_retention: {:?},
            is_intron_retained_in_frame: {:?},
            retains_rt_intron: {:?},
            has_rt_intron: {:?},
            has_rt_intron_map: {:?},
            retention_acceptor_score: {:?},
            retention_donor_score: {:?},
            ref_introns_component_size: {:?},
            query_component_size: {:?},
            component_retention_ratio: {:?},
            is_dirty_component: {:?},
            exonic_status: {:?},
            intronic_status: {:?}
            }}",
            self.intron_retention,
            self.retention_support_type,
            self.number_of_retentions,
            self.coords_of_retention,
            self.location_of_retention,
            self.is_intron_retained_in_frame,
            self.retains_rt_intron,
            self.has_rt_intron,
            self.has_rt_intron_map,
            self.retention_acceptor_score,
            self.retention_donor_score,
            self.ref_introns_component_size,
            self.query_component_size,
            self.component_retention_ratio,
            self.is_dirty_component,
            self.exonic_status,
            self.intronic_status
        )
    }
}

/// StartTruncationDescriptor struct
///
/// This struct represents the descriptor for
/// start truncation detection.
///
/// # Examples
///
/// ```rust, no_run
/// use isotools::modules::{StartTruncationDescriptor, StartTruncationValue};
///
/// let descriptor = StartTruncationDescriptor::new();
/// assert_eq!(descriptor.is_read_truncated, Value::Null);
/// ```
pub struct StartTruncationDescriptor {
    pub is_read_truncated: Value,
    pub is_novel_start: Value,
    pub is_dirty_component: Value,
    pub component_size: Value,
    pub ref_component_size: Value,
    pub query_component_size: Value,
    pub truncation_support_ratio: Value,
    pub is_truncation_supported: Value,
    pub component_truncation_ratio: Value,
}

/// StartTruncationDescriptor implementation
///
/// This implementation provides methods to create
/// a new instance of the StartTruncationDescriptor
impl StartTruncationDescriptor {
    /// Creates a new instance of the StartTruncationDescriptor.
    ///
    /// # Returns
    ///
    /// * A boxed instance of the StartTruncationDescriptor.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::StartTruncationDescriptor;
    ///
    /// let descriptor = StartTruncationDescriptor::new();
    /// assert_eq!(descriptor.is_read_truncated, Value::Null);
    /// ```
    pub fn new() -> Box<Self> {
        Box::new(Self {
            is_read_truncated: Value::Null,
            is_novel_start: Value::Null,
            is_dirty_component: Value::Null,
            component_size: Value::Null,
            ref_component_size: Value::Null,
            query_component_size: Value::Null,
            truncation_support_ratio: Value::Null,
            is_truncation_supported: Value::Null,
            component_truncation_ratio: Value::Null,
        })
    }
}

/// StartTruncationValue enum
///
/// This enum defines the keys used to access
/// the values in the StartTruncationDescriptor.
///
/// # Examples
///
/// ```rust, no_run
/// use isotools::modules::{StartTruncationDescriptor, StartTruncationValue};
///
/// let descriptor = StartTruncationDescriptor::new();
/// assert_eq!(descriptor.get_value(Box::new(StartTruncationValue::IsReadTruncated)), Some(Value::Null));
/// ```
#[derive(Debug, Clone)]
pub enum StartTruncationValue {
    IsReadTruncated,
    IsNovelStart,
    TruncationSupportRatio,
    IsTruncationSupported,
    ComponentSize,
    RefComponentSize,
    QueryComponentSize,
    ComponentTruncationRatio,
    IsDirtyComponent,
}

/// ModuleMap trait implementation for StartTruncationDescriptor
///
/// This implementation provides methods to get and set
/// values associated with the start truncation descriptor.
impl ModuleMap for StartTruncationDescriptor {
    /// Gets the value associated with the specified key.
    ///
    /// # Arguments
    ///
    /// * `key` - The key to retrieve the value for.
    ///
    /// # Returns
    ///
    /// * An `Option` containing the value associated with the key,
    ///  or `None` if the key is not found.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::{ModuleMap, StartTruncationDescriptor, StartTruncationValue};
    ///
    /// let module = StartTruncationDescriptor::new();
    /// assert_eq!(module.get_value(Box::new(StartTruncationValue::IsReadTruncated)), Some(Value::Null));
    /// ```
    fn get_value(&self, key: Box<dyn Any>) -> Option<serde_json::Value> {
        if let Ok(key) = key.downcast::<StartTruncationValue>() {
            match *key {
                StartTruncationValue::IsReadTruncated => Some(self.is_read_truncated.clone()),
                StartTruncationValue::IsNovelStart => Some(self.is_novel_start.clone()),
                StartTruncationValue::TruncationSupportRatio => {
                    Some(self.truncation_support_ratio.clone())
                }
                StartTruncationValue::IsTruncationSupported => {
                    Some(self.is_truncation_supported.clone())
                }
                StartTruncationValue::ComponentSize => Some(self.component_size.clone()),
                StartTruncationValue::RefComponentSize => Some(self.ref_component_size.clone()),
                StartTruncationValue::QueryComponentSize => Some(self.query_component_size.clone()),
                StartTruncationValue::ComponentTruncationRatio => {
                    Some(self.component_truncation_ratio.clone())
                }
                StartTruncationValue::IsDirtyComponent => Some(self.is_dirty_component.clone()),
            }
        } else {
            None
        }
    }

    /// Sets the value associated with the specified key.
    ///
    /// # Arguments
    ///
    /// * `key` - The key to set the value for.
    /// * `value` - The value to set.
    ///
    /// # Returns
    ///
    /// * A `Result` indicating success or failure.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::{ModuleMap, StartTruncationDescriptor, StartTruncationValue};
    ///
    /// let mut module = StartTruncationDescriptor::new();
    /// module.set_value(Box::new(StartTruncationValue::IsReadTruncated), Value::Bool(true)).unwrap();
    /// assert_eq!(module.get_value(Box::new(StartTruncationValue::IsReadTruncated)), Some(Value::Bool(true)));
    /// ```
    #[inline(always)]
    fn set_value(&mut self, key: Box<dyn Any>, value: Value) -> Result<(), String> {
        if let Ok(key) = key.downcast::<StartTruncationValue>() {
            match *key {
                StartTruncationValue::IsReadTruncated => {
                    self.is_read_truncated = value;
                    Ok(())
                }
                StartTruncationValue::IsNovelStart => {
                    self.is_novel_start = value;
                    Ok(())
                }
                StartTruncationValue::TruncationSupportRatio => {
                    self.truncation_support_ratio = value;
                    Ok(())
                }
                StartTruncationValue::IsTruncationSupported => {
                    self.is_truncation_supported = value;
                    Ok(())
                }
                StartTruncationValue::ComponentSize => {
                    self.component_size = value;
                    Ok(())
                }
                StartTruncationValue::RefComponentSize => {
                    self.ref_component_size = value;
                    Ok(())
                }
                StartTruncationValue::QueryComponentSize => {
                    self.query_component_size = value;
                    Ok(())
                }
                StartTruncationValue::ComponentTruncationRatio => {
                    self.component_truncation_ratio = value;
                    Ok(())
                }
                StartTruncationValue::IsDirtyComponent => {
                    self.is_dirty_component = value;
                    Ok(())
                }
            }
        } else {
            let err = format!("ERROR: You have tried to set a value for an unknown key!");
            log::error!("{}", err);
            Err(err)
        }
    }

    /// Downcasts the module descriptor to a specific type.
    ///
    /// # Returns
    ///
    /// * A reference to the module descriptor as a trait object.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::{ModuleMap, StartTruncationDescriptor};
    ///
    /// let module = StartTruncationDescriptor::new();
    /// assert_eq!(module.as_any().is::<StartTruncationDescriptor>(), true);
    /// ```
    fn as_any(&self) -> &dyn Any {
        self
    }

    /// Converts the module descriptor to JSON format.
    ///
    /// # Returns
    ///
    /// * A JSON object representing the module descriptor.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::{ModuleMap, StartTruncationDescriptor};
    ///
    /// let module = StartTruncationDescriptor::new();
    /// assert_eq!(module.to_json(), serde_json::json!({}));
    /// ```
    fn to_json(&self) -> Value {
        let mut map = Map::new();

        macro_rules! insert {
            ($field:ident) => {
                map.insert(stringify!($field).to_string(), self.$field.clone());
            };
        }

        insert!(is_read_truncated);
        insert!(is_novel_start);
        insert!(is_dirty_component);
        insert!(component_size);
        insert!(ref_component_size);
        insert!(query_component_size);
        insert!(truncation_support_ratio);
        insert!(is_truncation_supported);
        insert!(component_truncation_ratio);

        Value::Object(map)
    }
}

/// Debug implementation for StartTruncationDescriptor
///
/// This implementation provides a way to format the
/// StartTruncationDescriptor for debugging purposes.
impl std::fmt::Debug for StartTruncationDescriptor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{{
            is_read_truncated: {:?},
            is_novel_start: {:?},
            is_dirty_component: {:?},
            component_size: {:?},
            ref_component_size: {:?},
            query_component_size: {:?},
            truncation_support_ratio: {:?},
            is_truncation_supported: {:?},
            component_truncation_ratio: {:?}
            }}",
            self.is_read_truncated,
            self.is_novel_start,
            self.is_dirty_component,
            self.component_size,
            self.ref_component_size,
            self.query_component_size,
            self.truncation_support_ratio,
            self.is_truncation_supported,
            self.component_truncation_ratio
        )
    }
}

/// FusionDetectionDescriptor struct
///
/// This struct represents the descriptor for
/// fusion detection.
pub struct FusionDetectionDescriptor {
    is_fused_read: Value,
    is_fusion_supported: Value,
    component_size: Value,
    ref_component_size: Value,
    query_component_size: Value,
    whole_component_fusion_ratio: Value,
    real_component_fusion_ratio: Value,
    fake_component_fusion_ratio: Value,
    is_dirty_component: Value,
    location_of_fusion: Value,
    fusion_in_frame: Value,
}

/// FusionDetectionDescriptor implementation
///
/// This implementation provides methods to create
/// a new instance of the FusionDetectionDescriptor
impl FusionDetectionDescriptor {
    /// Creates a new instance of the FusionDetectionDescriptor.
    ///
    /// # Returns
    ///
    /// * A boxed instance of the FusionDetectionDescriptor.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::FusionDetectionDescriptor;
    ///
    /// let descriptor = FusionDetectionDescriptor::new();
    /// assert_eq!(descriptor.is_fused_read, Value::Null);
    /// ```
    pub fn new() -> Box<Self> {
        Box::new(Self {
            is_fused_read: Value::Null,
            is_fusion_supported: Value::Null,
            component_size: Value::Null,
            ref_component_size: Value::Null,
            query_component_size: Value::Null,
            whole_component_fusion_ratio: Value::Null,
            real_component_fusion_ratio: Value::Null,
            fake_component_fusion_ratio: Value::Null,
            is_dirty_component: Value::Null,
            location_of_fusion: Value::Null,
            fusion_in_frame: Value::Null,
        })
    }
}

/// FusionDetectionValue enum
///
/// This enum defines the keys used to access
/// the values in the FusionDetectionDescriptor.
pub enum FusionDetectionValue {
    IsFusedRead,
    IsFusionSupported,
    ComponentSize,
    RefComponentSize,
    QueryComponentSize,
    WholeComponentFusionRatio,
    RealComponentFusionRatio,
    FakeComponentFusionRatio,
    IsDirtyComponent,
    LocationOfFusion,
    FusionInFrame,
}

/// ModuleMap trait implementation for FusionDetectionDescriptor
impl ModuleMap for FusionDetectionDescriptor {
    /// Gets the value associated with the specified key.
    ///
    /// # Arguments
    ///
    /// * `key` - The key to retrieve the value for.
    ///
    /// # Returns
    ///
    /// * An `Option` containing the value associated with the key,
    ///  or `None` if the key is not found.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::{ModuleMap, FusionDetectionDescriptor, FusionDetectionValue};
    ///
    /// let module = FusionDetectionDescriptor::new();
    /// assert_eq!(module.get_value(Box::new(FusionDetectionValue::IsFusedRead)), Some(Value::Null));
    /// ```
    fn get_value(&self, key: Box<dyn Any>) -> Option<serde_json::Value> {
        if let Ok(key) = key.downcast::<FusionDetectionValue>() {
            match *key {
                FusionDetectionValue::IsFusedRead => Some(self.is_fused_read.clone()),
                FusionDetectionValue::IsFusionSupported => Some(self.is_fusion_supported.clone()),
                FusionDetectionValue::ComponentSize => Some(self.component_size.clone()),
                FusionDetectionValue::RefComponentSize => Some(self.ref_component_size.clone()),
                FusionDetectionValue::QueryComponentSize => Some(self.query_component_size.clone()),
                FusionDetectionValue::WholeComponentFusionRatio => {
                    Some(self.whole_component_fusion_ratio.clone())
                }
                FusionDetectionValue::RealComponentFusionRatio => {
                    Some(self.real_component_fusion_ratio.clone())
                }
                FusionDetectionValue::FakeComponentFusionRatio => {
                    Some(self.fake_component_fusion_ratio.clone())
                }
                FusionDetectionValue::IsDirtyComponent => Some(self.is_dirty_component.clone()),
                FusionDetectionValue::LocationOfFusion => Some(self.location_of_fusion.clone()),
                FusionDetectionValue::FusionInFrame => Some(self.fusion_in_frame.clone()),
            }
        } else {
            None
        }
    }

    /// Sets the value associated with the specified key.
    ///
    /// # Arguments
    ///
    /// * `key` - The key to set the value for.
    /// * `value` - The value to set.
    ///
    /// # Returns
    ///
    /// * A `Result` indicating success or failure.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::{ModuleMap, FusionDetectionDescriptor, FusionDetectionValue};
    ///
    /// let mut module = FusionDetectionDescriptor::new();
    /// module.set_value(Box::new(FusionDetectionValue::IsFusedRead), Value::Bool(true)).unwrap();
    /// assert_eq!(module.get_value(Box::new(FusionDetectionValue::IsFusedRead)), Some(Value::Bool(true)));
    /// ```
    #[inline(always)]
    fn set_value(&mut self, key: Box<dyn Any>, value: Value) -> Result<(), String> {
        if let Ok(key) = key.downcast::<FusionDetectionValue>() {
            match *key {
                FusionDetectionValue::IsFusedRead => {
                    self.is_fused_read = value;
                    Ok(())
                }
                FusionDetectionValue::IsFusionSupported => {
                    self.is_fusion_supported = value;
                    Ok(())
                }
                FusionDetectionValue::ComponentSize => {
                    self.component_size = value;
                    Ok(())
                }
                FusionDetectionValue::RefComponentSize => {
                    self.ref_component_size = value;
                    Ok(())
                }
                FusionDetectionValue::QueryComponentSize => {
                    self.query_component_size = value;
                    Ok(())
                }
                FusionDetectionValue::WholeComponentFusionRatio => {
                    self.whole_component_fusion_ratio = value;
                    Ok(())
                }
                FusionDetectionValue::RealComponentFusionRatio => {
                    self.real_component_fusion_ratio = value;
                    Ok(())
                }
                FusionDetectionValue::FakeComponentFusionRatio => {
                    self.fake_component_fusion_ratio = value;
                    Ok(())
                }
                FusionDetectionValue::IsDirtyComponent => {
                    self.is_dirty_component = value;
                    Ok(())
                }
                FusionDetectionValue::LocationOfFusion => {
                    self.location_of_fusion = value;
                    Ok(())
                }
                FusionDetectionValue::FusionInFrame => {
                    self.fusion_in_frame = value;
                    Ok(())
                }
            }
        } else {
            let err = format!("ERROR: You have tried to set a value for an unknown key!");
            log::error!("{}", err);
            Err(err)
        }
    }

    /// Downcasts the module descriptor to a specific type.
    ///
    /// # Returns
    ///
    /// * A reference to the module descriptor as a trait object.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::{ModuleMap, FusionDetectionDescriptor};
    ///
    /// let module = FusionDetectionDescriptor::new();
    /// assert_eq!(module.as_any().is::<FusionDetectionDescriptor>(), true);
    /// ```
    fn as_any(&self) -> &dyn Any {
        self
    }

    /// Converts the FusionDetectionDescriptor to a JSON object.
    ///
    /// # Returns
    ///
    /// * A JSON object representing the FusionDetectionDescriptor.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::FusionDetectionDescriptor;
    ///
    /// let descriptor = FusionDetectionDescriptor::new();
    /// let json = descriptor.to_json();
    /// assert_eq!(json["is_fused_read"], Value::Null);
    /// ```
    fn to_json(&self) -> Value {
        let mut map = Map::new();

        macro_rules! insert {
            ($field:ident) => {
                map.insert(stringify!($field).to_string(), self.$field.clone());
            };
        }

        insert!(is_fused_read);
        insert!(is_fusion_supported);
        insert!(component_size);
        insert!(ref_component_size);
        insert!(query_component_size);
        insert!(whole_component_fusion_ratio);
        insert!(real_component_fusion_ratio);
        insert!(fake_component_fusion_ratio);
        insert!(is_dirty_component);
        insert!(location_of_fusion);
        insert!(fusion_in_frame);

        Value::Object(map)
    }
}

/// Debug implementation for FusionDetectionDescriptor
///
/// This implementation provides a way to format the
/// FusionDetectionDescriptor for debugging purposes.
impl std::fmt::Debug for FusionDetectionDescriptor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{{
            is_fused_read: {:?},
            is_fusion_supported: {:?},
            component_size: {:?},
            ref_component_size: {:?},
            query_component_size: {:?}
            whole_component_fusion_ratio: {:?},
            real_component_fusion_ratio: {:?},
            fake_component_fusion_ratio: {:?},
            is_dirty_component: {:?},
            location_of_fusion: {:?},
            fusion_in_frame: {:?}
            }}",
            self.is_fused_read,
            self.is_fusion_supported,
            self.component_size,
            self.ref_component_size,
            self.query_component_size,
            self.whole_component_fusion_ratio,
            self.real_component_fusion_ratio,
            self.fake_component_fusion_ratio,
            self.is_dirty_component,
            self.location_of_fusion,
            self.fusion_in_frame,
        )
    }
}

/// PolyAPredictionDescriptor struct
///
/// This struct represents the descriptor for
/// poly A prediction detection.
pub struct PolyAPredictionDescriptor {
    pub is_poly_a_supported: Value,
    pub poly_a_score: Value,
    pub whole_poly_a_length: Value,
    pub genomic_poly_a: Value,
    pub is_intrapriming: Value,
    pub poly_a_location: Value,
    pub is_dirty_component: Value,
    pub intrapriming_comp_ratio: Value,
    pub forced_poly_a: Value,
}

/// PolyAPredictionDescriptor implementation
///
/// This implementation provides methods to create
/// a new instance of the PolyAPredictionDescriptor
impl PolyAPredictionDescriptor {
    /// Creates a new instance of the PolyAPredictionDescriptor.
    ///
    /// # Returns
    ///
    /// * A boxed instance of the PolyAPredictionDescriptor.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::PolyAPredictionDescriptor;
    ///
    /// let descriptor = PolyAPredictionDescriptor::new();
    /// assert_eq!(descriptor.is_poly_a, Value::Null);
    /// ```
    pub fn new() -> Box<Self> {
        Box::new(Self {
            is_poly_a_supported: Value::Null,
            poly_a_score: Value::Null,
            whole_poly_a_length: Value::Null,
            genomic_poly_a: Value::Null,
            is_intrapriming: Value::Null,
            poly_a_location: Value::Null,
            is_dirty_component: Value::Null,
            intrapriming_comp_ratio: Value::Null,
            forced_poly_a: Value::Null,
        })
    }
}

/// PolyAPredictionValue enum
///
/// This enum defines the keys used to access
/// the values in the PolyAPredictionDescriptor.
pub enum PolyAPredictionValue {
    IsPolyASupported,
    GenomicPolyA,
    PolyAScore,
    WholePolyALength,
    IsIntrapriming,
    PolyALocation,
    IsDirtyComponent,
    IntraprimingComponentRatio,
    ForcedPolyAPass,
}

/// ModuleMap trait implementation for PolyAPredictionDescriptor
///
/// This implementation provides methods to get and set
/// values associated with the poly A prediction descriptor.
impl ModuleMap for PolyAPredictionDescriptor {
    /// Gets the value associated with the specified key.
    ///
    /// # Arguments
    ///
    /// * `key` - The key to retrieve the value for.
    ///
    /// # Returns
    ///
    /// * An `Option` containing the value associated with the key,
    ///  or `None` if the key is not found.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::{ModuleMap, PolyAPredictionDescriptor, PolyAPredictionValue};
    ///
    /// let module = PolyAPredictionDescriptor::new();
    /// assert_eq!(module.get_value(Box::new(PolyAPredictionValue::IsPolyA)), Some(Value::Null));
    /// ```
    fn get_value(&self, key: Box<dyn Any>) -> Option<serde_json::Value> {
        if let Ok(key) = key.downcast::<PolyAPredictionValue>() {
            match *key {
                PolyAPredictionValue::IsPolyASupported => Some(self.is_poly_a_supported.clone()),
                PolyAPredictionValue::PolyAScore => Some(self.poly_a_score.clone()),
                PolyAPredictionValue::WholePolyALength => Some(self.whole_poly_a_length.clone()),
                PolyAPredictionValue::PolyALocation => Some(self.poly_a_location.clone()),
                PolyAPredictionValue::IsDirtyComponent => Some(self.is_dirty_component.clone()),
                PolyAPredictionValue::GenomicPolyA => Some(self.genomic_poly_a.clone()),
                PolyAPredictionValue::IsIntrapriming => Some(self.is_intrapriming.clone()),
                PolyAPredictionValue::IntraprimingComponentRatio => {
                    Some(self.intrapriming_comp_ratio.clone())
                }
                PolyAPredictionValue::ForcedPolyAPass => Some(Value::Null),
            }
        } else {
            None
        }
    }

    /// Sets the value associated with the specified key.
    ///
    /// # Arguments
    ///
    /// * `key` - The key to set the value for.
    /// * `value` - The value to set.
    ///
    /// # Returns
    ///
    /// * A `Result` indicating success or failure.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::{ModuleMap, PolyAPredictionDescriptor, PolyAPredictionValue};
    ///
    /// let mut module = PolyAPredictionDescriptor::new();
    /// module.set_value(Box::new(PolyAPredictionValue::IsPolyA), Value::Bool(true)).unwrap();
    /// assert_eq!(module.get_value(Box::new(PolyAPredictionValue::IsPolyA)), Some(Value::Bool(true)));
    /// ```
    #[inline(always)]
    fn set_value(&mut self, key: Box<dyn Any>, value: Value) -> Result<(), String> {
        if let Ok(key) = key.downcast::<PolyAPredictionValue>() {
            match *key {
                PolyAPredictionValue::IsPolyASupported => {
                    self.is_poly_a_supported = value;
                    Ok(())
                }
                PolyAPredictionValue::PolyAScore => {
                    self.poly_a_score = value;
                    Ok(())
                }
                PolyAPredictionValue::PolyALocation => {
                    self.poly_a_location = value;
                    Ok(())
                }
                PolyAPredictionValue::IsDirtyComponent => {
                    self.is_dirty_component = value;
                    Ok(())
                }
                PolyAPredictionValue::WholePolyALength => {
                    self.whole_poly_a_length = value;
                    Ok(())
                }
                PolyAPredictionValue::GenomicPolyA => {
                    self.genomic_poly_a = value;
                    Ok(())
                }
                PolyAPredictionValue::IsIntrapriming => {
                    self.is_intrapriming = value;
                    Ok(())
                }
                PolyAPredictionValue::IntraprimingComponentRatio => {
                    self.intrapriming_comp_ratio = value;
                    Ok(())
                }
                PolyAPredictionValue::ForcedPolyAPass => {
                    self.forced_poly_a = value;
                    Ok(())
                }
            }
        } else {
            let err = format!("ERROR: You have tried to set a value for an unknown key!");
            log::error!("{}", err);
            Err(err)
        }
    }

    /// Downcasts the module descriptor to a specific type.
    ///
    /// # Returns
    ///
    /// * A reference to the module descriptor as a trait object.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::{ModuleMap, PolyAPredictionDescriptor};
    ///
    /// let module = PolyAPredictionDescriptor::new();
    /// assert_eq!(module.as_any().is::<PolyAPredictionDescriptor>(), true);
    /// ```
    fn as_any(&self) -> &dyn Any {
        self
    }

    /// Converts the PolyAPredictionDescriptor to a JSON object.
    ///
    /// # Returns
    ///
    /// * A JSON object representing the PolyAPredictionDescriptor.
    ///
    /// # Examples
    ///
    /// ```rust, no_run
    /// use isotools::modules::PolyAPredictionDescriptor;
    ///
    /// let descriptor = PolyAPredictionDescriptor::new();
    /// let json = descriptor.to_json();
    /// assert_eq!(json["is_poly_a_supported"], Value::Null);
    /// ```
    fn to_json(&self) -> Value {
        let mut map = Map::new();

        macro_rules! insert {
            ($field:ident) => {
                map.insert(stringify!($field).to_string(), self.$field.clone());
            };
        }

        insert!(is_poly_a_supported);
        insert!(poly_a_score);
        insert!(whole_poly_a_length);
        insert!(genomic_poly_a);
        insert!(is_intrapriming);
        insert!(poly_a_location);
        insert!(is_dirty_component);
        insert!(intrapriming_comp_ratio);
        insert!(forced_poly_a);

        Value::Object(map)
    }
}

/// Debug implementation for PolyAPredictionDescriptor
///
/// This implementation provides a way to format the
/// PolyAPredictionDescriptor for debugging purposes.
impl std::fmt::Debug for PolyAPredictionDescriptor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{{
            is_poly_a_supported: {:?},
            poly_a_score: {:?},
            whole_poly_a_length: {:?},
            genomic_poly_a: {:?},
            is_intrapriming: {:?},
            poly_a_location: {:?},
            is_dirty_component: {:?},
            intrapriming_comp_ratio: {:?}
            forced_poly_a: {:?}
            }}",
            self.is_poly_a_supported,
            self.poly_a_score,
            self.whole_poly_a_length,
            self.genomic_poly_a,
            self.is_intrapriming,
            self.poly_a_location,
            self.is_dirty_component,
            self.intrapriming_comp_ratio,
            self.forced_poly_a
        )
    }
}
