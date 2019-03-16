// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: genomestore.proto

#ifndef PROTOBUF_INCLUDED_genomestore_2eproto
#define PROTOBUF_INCLUDED_genomestore_2eproto

#include <string>

#include <google/protobuf/stubs/common.h>

#if GOOGLE_PROTOBUF_VERSION < 3006001
#error This file was generated by a newer version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please update
#error your headers.
#endif
#if 3006001 < GOOGLE_PROTOBUF_MIN_PROTOC_VERSION
#error This file was generated by an older version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please
#error regenerate this file with a newer version of protoc.
#endif

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/arena.h>
#include <google/protobuf/arenastring.h>
#include <google/protobuf/generated_message_table_driven.h>
#include <google/protobuf/generated_message_util.h>
#include <google/protobuf/inlined_string_field.h>
#include <google/protobuf/metadata.h>
#include <google/protobuf/message.h>
#include <google/protobuf/repeated_field.h>  // IWYU pragma: export
#include <google/protobuf/extension_set.h>  // IWYU pragma: export
#include <google/protobuf/generated_enum_reflection.h>
#include <google/protobuf/unknown_field_set.h>
// @@protoc_insertion_point(includes)
#define PROTOBUF_INTERNAL_EXPORT_protobuf_genomestore_2eproto 

namespace protobuf_genomestore_2eproto {
// Internal implementation detail -- do not use these members.
struct TableStruct {
  static const ::google::protobuf::internal::ParseTableField entries[];
  static const ::google::protobuf::internal::AuxillaryParseTableField aux[];
  static const ::google::protobuf::internal::ParseTable schema[2];
  static const ::google::protobuf::internal::FieldMetadata field_metadata[];
  static const ::google::protobuf::internal::SerializationTable serialization_table[];
  static const ::google::protobuf::uint32 offsets[];
};
void AddDescriptors();
}  // namespace protobuf_genomestore_2eproto
namespace genomestore {
class Gene;
class GeneDefaultTypeInternal;
extern GeneDefaultTypeInternal _Gene_default_instance_;
class Interval;
class IntervalDefaultTypeInternal;
extern IntervalDefaultTypeInternal _Interval_default_instance_;
}  // namespace genomestore
namespace google {
namespace protobuf {
template<> ::genomestore::Gene* Arena::CreateMaybeMessage<::genomestore::Gene>(Arena*);
template<> ::genomestore::Interval* Arena::CreateMaybeMessage<::genomestore::Interval>(Arena*);
}  // namespace protobuf
}  // namespace google
namespace genomestore {

enum Base {
  a = 0,
  c = 1,
  g = 2,
  t = 3,
  u = 4,
  n = 5,
  Base_INT_MIN_SENTINEL_DO_NOT_USE_ = ::google::protobuf::kint32min,
  Base_INT_MAX_SENTINEL_DO_NOT_USE_ = ::google::protobuf::kint32max
};
bool Base_IsValid(int value);
const Base Base_MIN = a;
const Base Base_MAX = n;
const int Base_ARRAYSIZE = Base_MAX + 1;

const ::google::protobuf::EnumDescriptor* Base_descriptor();
inline const ::std::string& Base_Name(Base value) {
  return ::google::protobuf::internal::NameOfEnum(
    Base_descriptor(), value);
}
inline bool Base_Parse(
    const ::std::string& name, Base* value) {
  return ::google::protobuf::internal::ParseNamedEnum<Base>(
    Base_descriptor(), name, value);
}
// ===================================================================

class Interval : public ::google::protobuf::Message /* @@protoc_insertion_point(class_definition:genomestore.Interval) */ {
 public:
  Interval();
  virtual ~Interval();

  Interval(const Interval& from);

  inline Interval& operator=(const Interval& from) {
    CopyFrom(from);
    return *this;
  }
  #if LANG_CXX11
  Interval(Interval&& from) noexcept
    : Interval() {
    *this = ::std::move(from);
  }

  inline Interval& operator=(Interval&& from) noexcept {
    if (GetArenaNoVirtual() == from.GetArenaNoVirtual()) {
      if (this != &from) InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }
  #endif
  static const ::google::protobuf::Descriptor* descriptor();
  static const Interval& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const Interval* internal_default_instance() {
    return reinterpret_cast<const Interval*>(
               &_Interval_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  void Swap(Interval* other);
  friend void swap(Interval& a, Interval& b) {
    a.Swap(&b);
  }

  // implements Message ----------------------------------------------

  inline Interval* New() const final {
    return CreateMaybeMessage<Interval>(NULL);
  }

  Interval* New(::google::protobuf::Arena* arena) const final {
    return CreateMaybeMessage<Interval>(arena);
  }
  void CopyFrom(const ::google::protobuf::Message& from) final;
  void MergeFrom(const ::google::protobuf::Message& from) final;
  void CopyFrom(const Interval& from);
  void MergeFrom(const Interval& from);
  void Clear() final;
  bool IsInitialized() const final;

  size_t ByteSizeLong() const final;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input) final;
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const final;
  ::google::protobuf::uint8* InternalSerializeWithCachedSizesToArray(
      bool deterministic, ::google::protobuf::uint8* target) const final;
  int GetCachedSize() const final { return _cached_size_.Get(); }

  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const final;
  void InternalSwap(Interval* other);
  private:
  inline ::google::protobuf::Arena* GetArenaNoVirtual() const {
    return NULL;
  }
  inline void* MaybeArenaPtr() const {
    return NULL;
  }
  public:

  ::google::protobuf::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  // repeated .genomestore.Base seq = 6;
  int seq_size() const;
  void clear_seq();
  static const int kSeqFieldNumber = 6;
  ::genomestore::Base seq(int index) const;
  void set_seq(int index, ::genomestore::Base value);
  void add_seq(::genomestore::Base value);
  const ::google::protobuf::RepeatedField<int>& seq() const;
  ::google::protobuf::RepeatedField<int>* mutable_seq();

  // string ref = 1;
  void clear_ref();
  static const int kRefFieldNumber = 1;
  const ::std::string& ref() const;
  void set_ref(const ::std::string& value);
  #if LANG_CXX11
  void set_ref(::std::string&& value);
  #endif
  void set_ref(const char* value);
  void set_ref(const char* value, size_t size);
  ::std::string* mutable_ref();
  ::std::string* release_ref();
  void set_allocated_ref(::std::string* ref);

  // string chr = 2;
  void clear_chr();
  static const int kChrFieldNumber = 2;
  const ::std::string& chr() const;
  void set_chr(const ::std::string& value);
  #if LANG_CXX11
  void set_chr(::std::string&& value);
  #endif
  void set_chr(const char* value);
  void set_chr(const char* value, size_t size);
  ::std::string* mutable_chr();
  ::std::string* release_chr();
  void set_allocated_chr(::std::string* chr);

  // string seqs = 7;
  void clear_seqs();
  static const int kSeqsFieldNumber = 7;
  const ::std::string& seqs() const;
  void set_seqs(const ::std::string& value);
  #if LANG_CXX11
  void set_seqs(::std::string&& value);
  #endif
  void set_seqs(const char* value);
  void set_seqs(const char* value, size_t size);
  ::std::string* mutable_seqs();
  ::std::string* release_seqs();
  void set_allocated_seqs(::std::string* seqs);

  // uint32 start = 3;
  void clear_start();
  static const int kStartFieldNumber = 3;
  ::google::protobuf::uint32 start() const;
  void set_start(::google::protobuf::uint32 value);

  // uint32 len = 4;
  void clear_len();
  static const int kLenFieldNumber = 4;
  ::google::protobuf::uint32 len() const;
  void set_len(::google::protobuf::uint32 value);

  // bool strand = 5;
  void clear_strand();
  static const int kStrandFieldNumber = 5;
  bool strand() const;
  void set_strand(bool value);

  // float score = 8;
  void clear_score();
  static const int kScoreFieldNumber = 8;
  float score() const;
  void set_score(float value);

  // @@protoc_insertion_point(class_scope:genomestore.Interval)
 private:

  ::google::protobuf::internal::InternalMetadataWithArena _internal_metadata_;
  ::google::protobuf::RepeatedField<int> seq_;
  mutable int _seq_cached_byte_size_;
  ::google::protobuf::internal::ArenaStringPtr ref_;
  ::google::protobuf::internal::ArenaStringPtr chr_;
  ::google::protobuf::internal::ArenaStringPtr seqs_;
  ::google::protobuf::uint32 start_;
  ::google::protobuf::uint32 len_;
  bool strand_;
  float score_;
  mutable ::google::protobuf::internal::CachedSize _cached_size_;
  friend struct ::protobuf_genomestore_2eproto::TableStruct;
};
// -------------------------------------------------------------------

class Gene : public ::google::protobuf::Message /* @@protoc_insertion_point(class_definition:genomestore.Gene) */ {
 public:
  Gene();
  virtual ~Gene();

  Gene(const Gene& from);

  inline Gene& operator=(const Gene& from) {
    CopyFrom(from);
    return *this;
  }
  #if LANG_CXX11
  Gene(Gene&& from) noexcept
    : Gene() {
    *this = ::std::move(from);
  }

  inline Gene& operator=(Gene&& from) noexcept {
    if (GetArenaNoVirtual() == from.GetArenaNoVirtual()) {
      if (this != &from) InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }
  #endif
  static const ::google::protobuf::Descriptor* descriptor();
  static const Gene& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const Gene* internal_default_instance() {
    return reinterpret_cast<const Gene*>(
               &_Gene_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    1;

  void Swap(Gene* other);
  friend void swap(Gene& a, Gene& b) {
    a.Swap(&b);
  }

  // implements Message ----------------------------------------------

  inline Gene* New() const final {
    return CreateMaybeMessage<Gene>(NULL);
  }

  Gene* New(::google::protobuf::Arena* arena) const final {
    return CreateMaybeMessage<Gene>(arena);
  }
  void CopyFrom(const ::google::protobuf::Message& from) final;
  void MergeFrom(const ::google::protobuf::Message& from) final;
  void CopyFrom(const Gene& from);
  void MergeFrom(const Gene& from);
  void Clear() final;
  bool IsInitialized() const final;

  size_t ByteSizeLong() const final;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input) final;
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const final;
  ::google::protobuf::uint8* InternalSerializeWithCachedSizesToArray(
      bool deterministic, ::google::protobuf::uint8* target) const final;
  int GetCachedSize() const final { return _cached_size_.Get(); }

  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const final;
  void InternalSwap(Gene* other);
  private:
  inline ::google::protobuf::Arena* GetArenaNoVirtual() const {
    return NULL;
  }
  inline void* MaybeArenaPtr() const {
    return NULL;
  }
  public:

  ::google::protobuf::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  // repeated .genomestore.Interval exons = 8;
  int exons_size() const;
  void clear_exons();
  static const int kExonsFieldNumber = 8;
  ::genomestore::Interval* mutable_exons(int index);
  ::google::protobuf::RepeatedPtrField< ::genomestore::Interval >*
      mutable_exons();
  const ::genomestore::Interval& exons(int index) const;
  ::genomestore::Interval* add_exons();
  const ::google::protobuf::RepeatedPtrField< ::genomestore::Interval >&
      exons() const;

  // repeated .genomestore.Interval introns = 9;
  int introns_size() const;
  void clear_introns();
  static const int kIntronsFieldNumber = 9;
  ::genomestore::Interval* mutable_introns(int index);
  ::google::protobuf::RepeatedPtrField< ::genomestore::Interval >*
      mutable_introns();
  const ::genomestore::Interval& introns(int index) const;
  ::genomestore::Interval* add_introns();
  const ::google::protobuf::RepeatedPtrField< ::genomestore::Interval >&
      introns() const;

  // string id = 1;
  void clear_id();
  static const int kIdFieldNumber = 1;
  const ::std::string& id() const;
  void set_id(const ::std::string& value);
  #if LANG_CXX11
  void set_id(::std::string&& value);
  #endif
  void set_id(const char* value);
  void set_id(const char* value, size_t size);
  ::std::string* mutable_id();
  ::std::string* release_id();
  void set_allocated_id(::std::string* id);

  // uint32 tx_start = 2;
  void clear_tx_start();
  static const int kTxStartFieldNumber = 2;
  ::google::protobuf::uint32 tx_start() const;
  void set_tx_start(::google::protobuf::uint32 value);

  // uint32 tx_end = 3;
  void clear_tx_end();
  static const int kTxEndFieldNumber = 3;
  ::google::protobuf::uint32 tx_end() const;
  void set_tx_end(::google::protobuf::uint32 value);

  // uint32 cds_start = 4;
  void clear_cds_start();
  static const int kCdsStartFieldNumber = 4;
  ::google::protobuf::uint32 cds_start() const;
  void set_cds_start(::google::protobuf::uint32 value);

  // uint32 cds_end = 5;
  void clear_cds_end();
  static const int kCdsEndFieldNumber = 5;
  ::google::protobuf::uint32 cds_end() const;
  void set_cds_end(::google::protobuf::uint32 value);

  // uint32 idx = 6;
  void clear_idx();
  static const int kIdxFieldNumber = 6;
  ::google::protobuf::uint32 idx() const;
  void set_idx(::google::protobuf::uint32 value);

  // bool strand = 7;
  void clear_strand();
  static const int kStrandFieldNumber = 7;
  bool strand() const;
  void set_strand(bool value);

  // @@protoc_insertion_point(class_scope:genomestore.Gene)
 private:

  ::google::protobuf::internal::InternalMetadataWithArena _internal_metadata_;
  ::google::protobuf::RepeatedPtrField< ::genomestore::Interval > exons_;
  ::google::protobuf::RepeatedPtrField< ::genomestore::Interval > introns_;
  ::google::protobuf::internal::ArenaStringPtr id_;
  ::google::protobuf::uint32 tx_start_;
  ::google::protobuf::uint32 tx_end_;
  ::google::protobuf::uint32 cds_start_;
  ::google::protobuf::uint32 cds_end_;
  ::google::protobuf::uint32 idx_;
  bool strand_;
  mutable ::google::protobuf::internal::CachedSize _cached_size_;
  friend struct ::protobuf_genomestore_2eproto::TableStruct;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// Interval

// string ref = 1;
inline void Interval::clear_ref() {
  ref_.ClearToEmptyNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline const ::std::string& Interval::ref() const {
  // @@protoc_insertion_point(field_get:genomestore.Interval.ref)
  return ref_.GetNoArena();
}
inline void Interval::set_ref(const ::std::string& value) {
  
  ref_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), value);
  // @@protoc_insertion_point(field_set:genomestore.Interval.ref)
}
#if LANG_CXX11
inline void Interval::set_ref(::std::string&& value) {
  
  ref_.SetNoArena(
    &::google::protobuf::internal::GetEmptyStringAlreadyInited(), ::std::move(value));
  // @@protoc_insertion_point(field_set_rvalue:genomestore.Interval.ref)
}
#endif
inline void Interval::set_ref(const char* value) {
  GOOGLE_DCHECK(value != NULL);
  
  ref_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), ::std::string(value));
  // @@protoc_insertion_point(field_set_char:genomestore.Interval.ref)
}
inline void Interval::set_ref(const char* value, size_t size) {
  
  ref_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(),
      ::std::string(reinterpret_cast<const char*>(value), size));
  // @@protoc_insertion_point(field_set_pointer:genomestore.Interval.ref)
}
inline ::std::string* Interval::mutable_ref() {
  
  // @@protoc_insertion_point(field_mutable:genomestore.Interval.ref)
  return ref_.MutableNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline ::std::string* Interval::release_ref() {
  // @@protoc_insertion_point(field_release:genomestore.Interval.ref)
  
  return ref_.ReleaseNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline void Interval::set_allocated_ref(::std::string* ref) {
  if (ref != NULL) {
    
  } else {
    
  }
  ref_.SetAllocatedNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), ref);
  // @@protoc_insertion_point(field_set_allocated:genomestore.Interval.ref)
}

// string chr = 2;
inline void Interval::clear_chr() {
  chr_.ClearToEmptyNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline const ::std::string& Interval::chr() const {
  // @@protoc_insertion_point(field_get:genomestore.Interval.chr)
  return chr_.GetNoArena();
}
inline void Interval::set_chr(const ::std::string& value) {
  
  chr_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), value);
  // @@protoc_insertion_point(field_set:genomestore.Interval.chr)
}
#if LANG_CXX11
inline void Interval::set_chr(::std::string&& value) {
  
  chr_.SetNoArena(
    &::google::protobuf::internal::GetEmptyStringAlreadyInited(), ::std::move(value));
  // @@protoc_insertion_point(field_set_rvalue:genomestore.Interval.chr)
}
#endif
inline void Interval::set_chr(const char* value) {
  GOOGLE_DCHECK(value != NULL);
  
  chr_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), ::std::string(value));
  // @@protoc_insertion_point(field_set_char:genomestore.Interval.chr)
}
inline void Interval::set_chr(const char* value, size_t size) {
  
  chr_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(),
      ::std::string(reinterpret_cast<const char*>(value), size));
  // @@protoc_insertion_point(field_set_pointer:genomestore.Interval.chr)
}
inline ::std::string* Interval::mutable_chr() {
  
  // @@protoc_insertion_point(field_mutable:genomestore.Interval.chr)
  return chr_.MutableNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline ::std::string* Interval::release_chr() {
  // @@protoc_insertion_point(field_release:genomestore.Interval.chr)
  
  return chr_.ReleaseNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline void Interval::set_allocated_chr(::std::string* chr) {
  if (chr != NULL) {
    
  } else {
    
  }
  chr_.SetAllocatedNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), chr);
  // @@protoc_insertion_point(field_set_allocated:genomestore.Interval.chr)
}

// uint32 start = 3;
inline void Interval::clear_start() {
  start_ = 0u;
}
inline ::google::protobuf::uint32 Interval::start() const {
  // @@protoc_insertion_point(field_get:genomestore.Interval.start)
  return start_;
}
inline void Interval::set_start(::google::protobuf::uint32 value) {
  
  start_ = value;
  // @@protoc_insertion_point(field_set:genomestore.Interval.start)
}

// uint32 len = 4;
inline void Interval::clear_len() {
  len_ = 0u;
}
inline ::google::protobuf::uint32 Interval::len() const {
  // @@protoc_insertion_point(field_get:genomestore.Interval.len)
  return len_;
}
inline void Interval::set_len(::google::protobuf::uint32 value) {
  
  len_ = value;
  // @@protoc_insertion_point(field_set:genomestore.Interval.len)
}

// bool strand = 5;
inline void Interval::clear_strand() {
  strand_ = false;
}
inline bool Interval::strand() const {
  // @@protoc_insertion_point(field_get:genomestore.Interval.strand)
  return strand_;
}
inline void Interval::set_strand(bool value) {
  
  strand_ = value;
  // @@protoc_insertion_point(field_set:genomestore.Interval.strand)
}

// repeated .genomestore.Base seq = 6;
inline int Interval::seq_size() const {
  return seq_.size();
}
inline void Interval::clear_seq() {
  seq_.Clear();
}
inline ::genomestore::Base Interval::seq(int index) const {
  // @@protoc_insertion_point(field_get:genomestore.Interval.seq)
  return static_cast< ::genomestore::Base >(seq_.Get(index));
}
inline void Interval::set_seq(int index, ::genomestore::Base value) {
  seq_.Set(index, value);
  // @@protoc_insertion_point(field_set:genomestore.Interval.seq)
}
inline void Interval::add_seq(::genomestore::Base value) {
  seq_.Add(value);
  // @@protoc_insertion_point(field_add:genomestore.Interval.seq)
}
inline const ::google::protobuf::RepeatedField<int>&
Interval::seq() const {
  // @@protoc_insertion_point(field_list:genomestore.Interval.seq)
  return seq_;
}
inline ::google::protobuf::RepeatedField<int>*
Interval::mutable_seq() {
  // @@protoc_insertion_point(field_mutable_list:genomestore.Interval.seq)
  return &seq_;
}

// string seqs = 7;
inline void Interval::clear_seqs() {
  seqs_.ClearToEmptyNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline const ::std::string& Interval::seqs() const {
  // @@protoc_insertion_point(field_get:genomestore.Interval.seqs)
  return seqs_.GetNoArena();
}
inline void Interval::set_seqs(const ::std::string& value) {
  
  seqs_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), value);
  // @@protoc_insertion_point(field_set:genomestore.Interval.seqs)
}
#if LANG_CXX11
inline void Interval::set_seqs(::std::string&& value) {
  
  seqs_.SetNoArena(
    &::google::protobuf::internal::GetEmptyStringAlreadyInited(), ::std::move(value));
  // @@protoc_insertion_point(field_set_rvalue:genomestore.Interval.seqs)
}
#endif
inline void Interval::set_seqs(const char* value) {
  GOOGLE_DCHECK(value != NULL);
  
  seqs_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), ::std::string(value));
  // @@protoc_insertion_point(field_set_char:genomestore.Interval.seqs)
}
inline void Interval::set_seqs(const char* value, size_t size) {
  
  seqs_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(),
      ::std::string(reinterpret_cast<const char*>(value), size));
  // @@protoc_insertion_point(field_set_pointer:genomestore.Interval.seqs)
}
inline ::std::string* Interval::mutable_seqs() {
  
  // @@protoc_insertion_point(field_mutable:genomestore.Interval.seqs)
  return seqs_.MutableNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline ::std::string* Interval::release_seqs() {
  // @@protoc_insertion_point(field_release:genomestore.Interval.seqs)
  
  return seqs_.ReleaseNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline void Interval::set_allocated_seqs(::std::string* seqs) {
  if (seqs != NULL) {
    
  } else {
    
  }
  seqs_.SetAllocatedNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), seqs);
  // @@protoc_insertion_point(field_set_allocated:genomestore.Interval.seqs)
}

// float score = 8;
inline void Interval::clear_score() {
  score_ = 0;
}
inline float Interval::score() const {
  // @@protoc_insertion_point(field_get:genomestore.Interval.score)
  return score_;
}
inline void Interval::set_score(float value) {
  
  score_ = value;
  // @@protoc_insertion_point(field_set:genomestore.Interval.score)
}

// -------------------------------------------------------------------

// Gene

// string id = 1;
inline void Gene::clear_id() {
  id_.ClearToEmptyNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline const ::std::string& Gene::id() const {
  // @@protoc_insertion_point(field_get:genomestore.Gene.id)
  return id_.GetNoArena();
}
inline void Gene::set_id(const ::std::string& value) {
  
  id_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), value);
  // @@protoc_insertion_point(field_set:genomestore.Gene.id)
}
#if LANG_CXX11
inline void Gene::set_id(::std::string&& value) {
  
  id_.SetNoArena(
    &::google::protobuf::internal::GetEmptyStringAlreadyInited(), ::std::move(value));
  // @@protoc_insertion_point(field_set_rvalue:genomestore.Gene.id)
}
#endif
inline void Gene::set_id(const char* value) {
  GOOGLE_DCHECK(value != NULL);
  
  id_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), ::std::string(value));
  // @@protoc_insertion_point(field_set_char:genomestore.Gene.id)
}
inline void Gene::set_id(const char* value, size_t size) {
  
  id_.SetNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(),
      ::std::string(reinterpret_cast<const char*>(value), size));
  // @@protoc_insertion_point(field_set_pointer:genomestore.Gene.id)
}
inline ::std::string* Gene::mutable_id() {
  
  // @@protoc_insertion_point(field_mutable:genomestore.Gene.id)
  return id_.MutableNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline ::std::string* Gene::release_id() {
  // @@protoc_insertion_point(field_release:genomestore.Gene.id)
  
  return id_.ReleaseNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited());
}
inline void Gene::set_allocated_id(::std::string* id) {
  if (id != NULL) {
    
  } else {
    
  }
  id_.SetAllocatedNoArena(&::google::protobuf::internal::GetEmptyStringAlreadyInited(), id);
  // @@protoc_insertion_point(field_set_allocated:genomestore.Gene.id)
}

// uint32 tx_start = 2;
inline void Gene::clear_tx_start() {
  tx_start_ = 0u;
}
inline ::google::protobuf::uint32 Gene::tx_start() const {
  // @@protoc_insertion_point(field_get:genomestore.Gene.tx_start)
  return tx_start_;
}
inline void Gene::set_tx_start(::google::protobuf::uint32 value) {
  
  tx_start_ = value;
  // @@protoc_insertion_point(field_set:genomestore.Gene.tx_start)
}

// uint32 tx_end = 3;
inline void Gene::clear_tx_end() {
  tx_end_ = 0u;
}
inline ::google::protobuf::uint32 Gene::tx_end() const {
  // @@protoc_insertion_point(field_get:genomestore.Gene.tx_end)
  return tx_end_;
}
inline void Gene::set_tx_end(::google::protobuf::uint32 value) {
  
  tx_end_ = value;
  // @@protoc_insertion_point(field_set:genomestore.Gene.tx_end)
}

// uint32 cds_start = 4;
inline void Gene::clear_cds_start() {
  cds_start_ = 0u;
}
inline ::google::protobuf::uint32 Gene::cds_start() const {
  // @@protoc_insertion_point(field_get:genomestore.Gene.cds_start)
  return cds_start_;
}
inline void Gene::set_cds_start(::google::protobuf::uint32 value) {
  
  cds_start_ = value;
  // @@protoc_insertion_point(field_set:genomestore.Gene.cds_start)
}

// uint32 cds_end = 5;
inline void Gene::clear_cds_end() {
  cds_end_ = 0u;
}
inline ::google::protobuf::uint32 Gene::cds_end() const {
  // @@protoc_insertion_point(field_get:genomestore.Gene.cds_end)
  return cds_end_;
}
inline void Gene::set_cds_end(::google::protobuf::uint32 value) {
  
  cds_end_ = value;
  // @@protoc_insertion_point(field_set:genomestore.Gene.cds_end)
}

// uint32 idx = 6;
inline void Gene::clear_idx() {
  idx_ = 0u;
}
inline ::google::protobuf::uint32 Gene::idx() const {
  // @@protoc_insertion_point(field_get:genomestore.Gene.idx)
  return idx_;
}
inline void Gene::set_idx(::google::protobuf::uint32 value) {
  
  idx_ = value;
  // @@protoc_insertion_point(field_set:genomestore.Gene.idx)
}

// bool strand = 7;
inline void Gene::clear_strand() {
  strand_ = false;
}
inline bool Gene::strand() const {
  // @@protoc_insertion_point(field_get:genomestore.Gene.strand)
  return strand_;
}
inline void Gene::set_strand(bool value) {
  
  strand_ = value;
  // @@protoc_insertion_point(field_set:genomestore.Gene.strand)
}

// repeated .genomestore.Interval exons = 8;
inline int Gene::exons_size() const {
  return exons_.size();
}
inline void Gene::clear_exons() {
  exons_.Clear();
}
inline ::genomestore::Interval* Gene::mutable_exons(int index) {
  // @@protoc_insertion_point(field_mutable:genomestore.Gene.exons)
  return exons_.Mutable(index);
}
inline ::google::protobuf::RepeatedPtrField< ::genomestore::Interval >*
Gene::mutable_exons() {
  // @@protoc_insertion_point(field_mutable_list:genomestore.Gene.exons)
  return &exons_;
}
inline const ::genomestore::Interval& Gene::exons(int index) const {
  // @@protoc_insertion_point(field_get:genomestore.Gene.exons)
  return exons_.Get(index);
}
inline ::genomestore::Interval* Gene::add_exons() {
  // @@protoc_insertion_point(field_add:genomestore.Gene.exons)
  return exons_.Add();
}
inline const ::google::protobuf::RepeatedPtrField< ::genomestore::Interval >&
Gene::exons() const {
  // @@protoc_insertion_point(field_list:genomestore.Gene.exons)
  return exons_;
}

// repeated .genomestore.Interval introns = 9;
inline int Gene::introns_size() const {
  return introns_.size();
}
inline void Gene::clear_introns() {
  introns_.Clear();
}
inline ::genomestore::Interval* Gene::mutable_introns(int index) {
  // @@protoc_insertion_point(field_mutable:genomestore.Gene.introns)
  return introns_.Mutable(index);
}
inline ::google::protobuf::RepeatedPtrField< ::genomestore::Interval >*
Gene::mutable_introns() {
  // @@protoc_insertion_point(field_mutable_list:genomestore.Gene.introns)
  return &introns_;
}
inline const ::genomestore::Interval& Gene::introns(int index) const {
  // @@protoc_insertion_point(field_get:genomestore.Gene.introns)
  return introns_.Get(index);
}
inline ::genomestore::Interval* Gene::add_introns() {
  // @@protoc_insertion_point(field_add:genomestore.Gene.introns)
  return introns_.Add();
}
inline const ::google::protobuf::RepeatedPtrField< ::genomestore::Interval >&
Gene::introns() const {
  // @@protoc_insertion_point(field_list:genomestore.Gene.introns)
  return introns_;
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__
// -------------------------------------------------------------------


// @@protoc_insertion_point(namespace_scope)

}  // namespace genomestore

namespace google {
namespace protobuf {

template <> struct is_proto_enum< ::genomestore::Base> : ::std::true_type {};
template <>
inline const EnumDescriptor* GetEnumDescriptor< ::genomestore::Base>() {
  return ::genomestore::Base_descriptor();
}

}  // namespace protobuf
}  // namespace google

// @@protoc_insertion_point(global_scope)

#endif  // PROTOBUF_INCLUDED_genomestore_2eproto